import json
from lenskappa.analysis import kappa_inference, kappaSet, inference
from pathlib import Path
import sys

def run_inference():
    config_path = sys.argv[1]
    with open(config_path) as f:
        config_data = json.load(f)
    
    base_inference = config_data["base-inference"]
    if base_inference == "kappa_set":
        mod = kappaSet
        path = Path(mod.__file__).parents[0] / "kappa_set_template.json"


    elif base_inference == "kappa":
        mod = kappa_inference
        path = Path(mod.__file__).parents[0] / "kappa_template.json"
    
    with open(path) as f:
        base_inference_config = json.load(f)
    
    inf = inference.build_inference(config_data, base_inference_config, mod)
    inf.run_inference()