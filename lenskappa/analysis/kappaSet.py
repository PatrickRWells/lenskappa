from pathlib import Path
from typing import List
from itertools import combinations
from lenskappa.analysis import kappa_inference
from lenskappa.analysis import inference
from lenskappa.analysis.transformation import Transformation
import json

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
        return analysis_object


class run_analyses(Transformation):
    def __call__(self, *args, **kwargs):
        return self.run_analyses(*args, **kwargs)
    def run_analyses(self, analyses: list):
        for analysis in analyses:
            analysis.run_inference()


