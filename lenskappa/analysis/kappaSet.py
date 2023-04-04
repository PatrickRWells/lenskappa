from pathlib import Path
from typing import List
from itertools import combinations
from lenskappa.analysis import kappa_inference
from lenskappa.analysis import inference
from lenskappa.analysis.transformation import Transformation
import json

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
    The list of weights to select from when building inferences.
nweights: int
    The number of weights from the weights list to use in each inference. For 
    each inference, the weights used will be any weights passed into "base-weights",
    plus nweights from the weights list. 
z_s: float
    The redshift of this lens. This is required for selecting the redshift plane
    to use for the weak lensing maps.

base_weights: list, default = None
    The list of weights that will be used for every inference in the set. Optional.


"""


class build_analyses(Transformation):
    def __call__(self, *args, **kwargs):
        return self.build_analyses(*args, **kwargs)
    def build_analyses(self,
        wnc_base_path: Path, ms_wnc_base_path: Path,
        wlm_base_path: Path, output_base_path: Path,
        wnc_base_names: List[str], weights: List[str], 
        nweights: int, z_s: float, base_weights: List[str] = None):
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
        for param_combo in weight_parameter_combinations:
            for combo in weight_combinations:
                analyses.append(self.build_single_analysis(*param_combo, wlm_base_path, combo, z_s))
        return analyses

    def build_single_analysis(self, wnc_path, ms_weight_paths, output_path, wlm_base_path, weight_combination, z_s):
        fname = "_".join(weight_combination) + ".k"
        new_output_path = Path(output_path) / fname


        parameters = {
            "base-inference": "kappa",
            "parameters": {
                "wnc_path": wnc_path,
                "ms_wnc_path": ms_weight_paths,
                "wlm_path": wlm_base_path,
                "weights": weight_combination,
                "z_s": z_s,
                "output_path": new_output_path
            }
        }
        kappa_module = kappa_inference
        mod_path = Path(kappa_module.__file__)
        template_path = mod_path.parents[0] / "kappa_template.json"
        with open(template_path, "r") as f:
            base_template = json.load(f)
        analysis_object = inference.build_inference(parameters, base_template , kappa_module)
        analysis_object.run_to("attach_ms_wlm")
        return analysis_object


class run_analyses(Transformation):
    def __call__(self, *args, **kwargs):
        return self.run_analyses(*args, **kwargs)
    def run_analyses(self, analyses: list):
        for analysis in analyses:
            analysis.run_inference()


