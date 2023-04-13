"""

Utilities for attaching weak lensing observables from maps to weighted number counts retrieved from
the milennium simulation

"""


from pathlib import Path
import heinlein
import lenskappa
import pandas as pd
from itertools import product
import numpy as np
import astropy.units as u
import multiprocessing as mp
from functools import partial
import re
import collections
import logging
from dask.distributed import get_client, secede, rejoin

logger = logging.getLogger("attach_wlm")


def attach_wlm(wnc: dict, redshift: float, redshift_plane: list = None, wlm_path: Path = None,  wnc_path: Path = None, inplace=True, threads = 1,*args, **kwargs):
    """
    Attaches weak lensing observables to a set of weighted number counts found
    in the millennium simulation. The kappa and gamma values are stored
    in a subfolder of the folder with the original values, making it easy 
    to re-attach them later.


    Parameters:
    -----------
    
    wnc_path: Path to the weighted number counts files. No particular filename is necessary, we just 
    expect somewhere in the name to find the tuple that tells us which field it is from
    
    redshift: Redshift cutoff used when computing weighted number counts

    wlm_path (optional): Path to the location of the weak lensing maps. If not provided,
    will look for them in heinlein

    inplace (optional): Whether or not to modify the weighte number count files in place. 
    If true, will add columns to the files. If false, will create new files.
    
    """
    if wlm_path is None:
        p = heinlein.get_path("ms", "wl_maps")
        if p is None:
            raise FileNotFoundError("No path found for weak lensing maps")
        wlm_path = p
    wnc_files = {}
    if wnc_path is not None:
        wnc_file_paths = [f for f in wnc_path.glob("*.csv")]
        file_fields = [re.search(r"\d_\d", f.name).group() for f in wnc_file_paths]
        file_fields = [tuple(int(val) for val in f.split("_")) for f in file_fields]
        wnc_files = {field: path for (field, path) in zip(file_fields, wnc_file_paths) if field in wnc.keys()}

        #To ensure order is the same for easier mapping
        wnc_files = {field: wnc_files[field] for field in wnc.keys()}

        


    if redshift_plane is None:
        plane_numbers = get_redshift_plane(redshift)
    elif type(redshift_plane) != list:
        plane_numbers = [redshift_plane]
    else:
        plane_numbers = redshift_plane
    plane_numbers.sort()
    f_ = partial(load_single_field, plane_numbers = plane_numbers, wlm_path = wlm_path)
    client = get_client()
    inputs = list(zip(wnc.values(), wnc_files.values(), wnc.keys()))
    
    map = client.map(f_, inputs)
    secede()
    wnc_dfs = client.gather(map)
    rejoin()
    if any([df is not None for df in wnc_dfs]):
        return pd.concat([df for df in wnc_dfs if df is not None])
    return pd.DataFrame()

def get_redshift_plane(redshift):
    """
    Figures out the best redshift plane to use given a particular redshift.
    This is only necesesary if the redshift plane is not passed into the 
    original attach_wlm function.
    """


    main_path = Path(lenskappa.__file__).parents[0]
    ms_distances_path = main_path / "datasets" / "surveys" / "ms" / "Millennium_distances.txt"
    ms_distances = pd.read_csv(ms_distances_path, delim_whitespace=True)
    ms_redshifts = ms_distances["Redshift"]
    closest_above = ms_redshifts[ms_redshifts > redshift].min()
    closest_below = ms_redshifts[ms_redshifts < redshift].max()

    if abs((redshift - closest_above) / redshift) <= 0.05:
        plane_redshift = closest_above

    elif abs((redshift - closest_below) / redshift) < 0.05:
        plane_redshift = closest_below

    while True:
        print(f"The closest redshift planes to z = {redshift} are z = {closest_above} and z = {closest_below}")
        print(f"Would you like to use one of these, or average the planes?")
        choice = input(f"Closest (A)bove ({closest_above}), closest (B)elow ({closest_below}), or avera(G)e: ").upper()
        if choice not in ["A", "B", "G"]:
            print("Invalid input")
            continue
        if choice == "A":
            plane_redshift = [closest_above]
        elif choice == "B":
            plane_redshift = [closest_below]
        elif choice == "G":
            plane_redshift = [closest_above, closest_below]
        break
    plane_numbers = ms_distances["PlaneNumber"][ms_redshifts.isin(plane_redshift)]
    return [pn for pn in plane_numbers]



def load_single_field(input_tuple, plane_numbers, wlm_path):
        """"
        The millennium simulation is split 64 4x4 deg^2 fields. It's assumed that the
        weighted number counts are similarly split up. This function loads a single 
        one of these fields for a given redshift plane (or planes, if using an average)
        
        """
        wnc = input_tuple[0]
        wnc_file = input_tuple[1]
        field = input_tuple[2]
        output_folder = wnc_file.parents[0] / "+".join([str(p) for p in plane_numbers])

        if check_for_wlm(wnc_file, plane_numbers):
            field_file = output_folder / wnc_file.name
            wlm_data = pd.read_csv(field_file)
            if len(wlm_data) == len(wnc):
                return pd.merge(wnc, wlm_data, on = ["ra", "dec"])



        output = wnc[["ra", "dec"]].copy()
        logger.info(f"Working on field {field}")
        wl_data = load_field_wlm(field, plane_numbers, wlm_path)
        if wl_data is None:
            return
        coords = list(zip(wnc["ra"], wnc["dec"]))
        indices = [get_index_from_position(*c) for c in coords]
        index_arrays = tuple(map(list, zip(*indices)))
        if not all([np.any(d["kappa"]) for d in wl_data.values()]):
            print(f"Warning: Tried to average two plains but they do not both have weak lensing maps for field {field}")
            return
        kappa_total = np.sum([d['kappa'] for d in wl_data.values()], axis = 0) / len(wl_data)
        gamma_total = np.sum([d['gamma'] for d in wl_data.values()], axis = 0) / len(wl_data)
        if not kappa_total.any():
            return
        kappas = kappa_total[index_arrays[0], index_arrays[1]]
        gammas = gamma_total[index_arrays[0], index_arrays[1]]

        output["kappa"] = kappas
        output["gamma"] = gammas
        output_folder.mkdir(exist_ok = True)
        output_file = output_folder / wnc_file.name
        output.to_csv(output_file, index=False)
        wnc["kappa"] = kappas
        wnc["gamma"] = gammas
        return wnc

def load_field_wlm(field_label, plane_numbers, wlm_path):
    """
    This function loads the weak lensing maps for a single field in a given redshift
    plane (or planes, if using an average). 
    
    """

    possible_filetypes = [".kappa", ".gamma_1", ".gamma_2", ".Phi"]
    searchname = "*N_4096_ang_4_rays_to_plane_29_f.*"
    files = [f for f in wlm_path.glob(searchname) if f.suffix in possible_filetypes]
    data = {}
    if not files:
        dirs = [path for path in wlm_path.glob("*") if path.is_dir()]
        dirnames = [path.name for path in dirs]
        if all(k in dirnames for k in ["kappa", "gamma"]):
            kappa = load_kappa(field_label, wlm_path / "kappa")
            gamma = load_gamma(field_label, wlm_path / "gamma")
            return {"kappa": kappa, "gamma": gamma}
        elif all([any([str(pn) in name for name in dirnames]) for pn in plane_numbers]):
            #Found a directory for all the requested planes
            if not check_for_wlm_maps(plane_numbers, field_label, wlm_path):
                return
            for pn in plane_numbers:
                plane_dir = [d for d in dirs if str(pn) in d.name]
                if len(plane_dir) != 1:
                    print(f"Found multiple directories for plane {pn}!")
                    exit()
                data.update({pn: load_field_wlm(field_label, [pn], plane_dir[0])})
            return data
        else:
            pns = ",".join([str(pn) for pn in plane_numbers])
            print(f"No .kappa and .gamma files found for field {field_label} in plane(s) {pns}")
            return

def check_for_wlm_maps(plane_numbers, field_label, wlm_path):
    """
    Checks to see if weak lensing maps exist for the given field for all
    redshift planes passed to this function.

    Returns True/False
    """
    dirs = [path for path in wlm_path.glob("*") if path.is_dir()]
    basename = f"GGL_los_8_{field_label[0]}_{field_label[1]}_N_4096_ang_4_rays_to_plane"
    for pn in plane_numbers:        
        plane_dir = [d for d in dirs if str(pn) in d.name]
        if len(plane_dir) != 1:
            print(f"Found multiple directories for plane {pn}!")
            return False
        plane_path = wlm_path / plane_dir[0]
        kappa_path = plane_path / "kappa"
        files = [file for file in kappa_path.glob("*.kappa") if basename in file.stem]
        if len(files) != 1:
            print(f"Found wrong number of files for kappa map for field {field_label}")
            return False
    
    return True

def check_for_wlm(wnc_path: Path, plane_numbers: list):
    """
    Checks to see if a given set of weighted number counts has
    already attached weak lensing values for the given plane numbers. 
    """
    plane_numbers.sort()
    fname = wnc_path.name
    folder = wnc_path.parents[0] / "+".join([str(p) for p in plane_numbers])
    file_path = folder / fname
    return file_path.exists()

def load_kappa(field_label, path):
    """
    Loads kappa files.
    """
    basename = f"GGL_los_8_{field_label[0]}_{field_label[1]}_N_4096_ang_4_rays_to_plane"
    files = [file for file in path.glob("*.kappa") if basename in file.stem]
    if len(files) != 1:
        print(f"Found wrong number of files for kappa map for field {field_label}")
        return []
    kf = files[0]
    kappa =  np.fromfile(kf, np.float32).reshape((4096, 4096))
    return kappa

def load_gamma(field_label, path):
    """
    Loads gamma files. x`
    """
    basename = f"GGL_los_8_{field_label[0]}_{field_label[1]}_N_4096_ang_4_rays_to_plane"

    gamma1_files = [file for file in path.glob("*.gamma_1") if basename in file.stem]
    gamma2_files = [file for file in path.glob("*.gamma_2") if basename in file.stem]
    if len(gamma1_files) != 1 or len(gamma2_files) != 1:
        print(f"Found wrong number of gamma files for field {field_label}")
        return []
    gamma1 = np.fromfile(gamma1_files[0], np.float32).reshape((4096, 4096))
    gamma2 = np.fromfile(gamma2_files[0], np.float32).reshape((4096, 4096))
    return np.sqrt(gamma1**2 + gamma2**2)


def get_index_from_position(pos_x, pos_y):
    """
    Returns the index of the nearest grid point given an angular position.
    
    """

    if pos_x > 2:
        pos_x = pos_x - 360
    pos_x = pos_x*u.deg
    pos_y = pos_y*u.deg

    l_field = (4.0*u.degree)
    n_pix = 4096.0
    l_pix = l_field/n_pix

    x_pix = (pos_x + 2.0*u.deg)/l_pix - 0.5
    y_pix = (pos_y + 2.0*u.deg)/l_pix - 0.5
    return int(round(x_pix.value)), int(round(y_pix.value))
