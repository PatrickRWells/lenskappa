from pathlib import Path
from typing import List
from itertools import combinations
from lenskappa.analysis import kappa_analysis
from lenskappa.analysis import analysis
from lenskappa.analysis.transformation import Transformation
from lenskappa.utils import attach_ms_wlm
import json
import logging
"""
A kappa set is a series of kappa inferences done on a single system.
The analysis can be done with several combinations of weights.

Parameters
----------
wnc_base_names: list
    The base names of the combination of limiting magnitude and aperture
    being considered. These names should be formatted as {limmag}_{aperture},
    and will be used to discover the relevant files for the field itself and the
    millennium simulation.
wnc_base_path: Path
    The base path where the weighted number counts for the system are located.
    They should be csvs with the same basenames discussed above. 
ms_wnc_base_path: Path
    The base path where the weighted number counts from the simulation are located.
    They should be placed in folders with the same base name discussed above
wlm_base_path: Path
    The base path where the weak lensing maps are located. Should be put in folders
    by redshift plane, named "PlaneXX" where XX is the plane number.
output_base_path: Path
    The location to place the outputs. Outputs will be placed in folders with the
    same base name discussed above, without file names dependent on the weights
    being used.
weights: list
    The list of weights to select from when building analysiss.
nweights: int
    The number of weights from the weights list to use in each analysis. For 
    each analysis, the weights used will be any weights passed into "base-weights",
    plus nweights from the weights list. 
z_s: float
    The redshift of this lens. This is required for selecting the redshift plane
    to use for the weak lensing maps.

base_weights: list, default = None
    The list of weights that will be used for every analysis in the set. Optional.


"""

class build_analyses(Transformation):
    def __call__(self, *args, **kwargs):
        return self.build_analyses(*args, **kwargs)
    
    def build_analyses(self,
        wnc_base_path: Path, ms_wnc_base_path: Path,
        wlm_base_path: Path, output_base_path: Path,
        wnc_base_names: List[str], weights: List[str], 
        nweights: int, z_s: float, logger_,  base_weights: List[str] = None, **kwargs):

        wnc_paths = [Path(wnc_base_path) / f"{bn}.csv" for bn in wnc_base_names]
        ms_weight_paths = [Path(ms_wnc_base_path) / f"{bn}" for bn in wnc_base_names]
        output_paths = [Path(output_base_path)  / bn for bn in wnc_base_names]
        for op in output_paths:
            op.mkdir(exist_ok=True, parents=True)

        weight_combinations = [list(c) for c in combinations(weights, nweights)]

        if base_weights is not None:
            if type(base_weights) != list:
                base_weights = [base_weights]
            weight_combinations = [base_weights + wc for wc in weight_combinations]
        weight_parameter_combinations = list(zip(wnc_paths, ms_weight_paths, output_paths))
        analyses = []
        self.plane = self.get_planes(z_s)
        for param_combo in weight_parameter_combinations:
            for combo in weight_combinations:
                analyses.append(self.build_single_analysis(*param_combo, wlm_base_path, combo, z_s, **kwargs))

        return analyses

    def build_single_analysis(self, wnc_path, ms_weight_paths, output_path, wlm_base_path, weight_combination, z_s, name, **kwargs):
        fname = "_".join(weight_combination) + ".k"
        new_output_path = Path(output_path) / fname


        parameters = {
            "base-analysis": "kappa",
            "name": name,
            "parameters": {
                "wnc_path": wnc_path,
                "ms_wnc_path": ms_weight_paths,
                "wlm_path": wlm_base_path,
                "weights": weight_combination,
                "z_s": z_s,
                "output_path": new_output_path,
                "redshift_plane": self.plane
            }
        }
        kappa_module = kappa_analysis
        mod_path = Path(kappa_module.__file__)
        template_path = mod_path.parents[0] / "kappa_template.json"
        with open(template_path, "r") as f:
            base_template = json.load(f)

        analysis_object = analysis.build_analysis(parameters, base_template , kappa_module, **kwargs)
        return analysis_object
    
    def get_planes(self, z_s: float):
        return attach_ms_wlm.get_redshift_plane(z_s)


