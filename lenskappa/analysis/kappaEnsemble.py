from pathlib import Path
from itertools import combinations
from lenskappa.analysis.transformation import Transformation
from lenskappa.analysis import kappaSet
from lenskappa.analysis.analysis import build_analysis
import json
import toml
from typing import Union
from copy import deepcopy
import logging
from loguru import logger
"""
A kappa ensemble is a kappa set that is being applied to
several different systems. When running a kappa anlysis, we
first have to go through the weighted number counts and extract
the weak lensing observables for each line of sight included
in those counts. This requires user input, so we go through
this process first. 


Parameters
----------

lens_parameters: dict
    A dictionary of lenses where the key is the name of the lens and the
    value is another dictionary with the following necessary paramters:
        wnc_base_path: The base path where the weighted number counts for this
            lens system are stored
        ms_wnc_base_path: The base path where the millennium simulation weighted number
            counts are stored for this lens
        z_s: The redshift of the source
    The remaining parameters will be applied to each lens, and should be provided in
    kappa_set_parameters (see below)
    
kappa_set_parameters: dict or Path
    The parameters for the kappa set to be run on each lens (minus the ones already included
    in the "lens_paramters" dictionary). This can be a dictionary with the actual parameter values, 
    or a path to a file that contains the dictionary (in json or toml format). See the documentation
    for the "kappa set" for details fo these parameters.

base_output_path: Path
    The location to put the outputs. All outputs will be placed in a folder with the same
    name as the lenses. Further subidivsion as discussed in the documentation for the "kappa-set"
    analysis. 
"""


class build_analyses(Transformation):
    """
    Builds the analyses
    
    """
    def __call__(self, logger_, *args, **kwargs):

        return self.build_analyses(logger_=logger_, *args, **kwargs)
    def build_analyses(self,
        lens_parameters: dict, kappa_set_parameters: Union[dict, Path],
        output_base_path: Path, logger_ = logging.Logger, **kwargs):
        if type(kappa_set_parameters) == dict:
            self.common_parameters = kappa_set_parameters
        else:
            with open(kappa_set_parameters, "r") as f:
                try:
                    self.common_parameters = json.load(f)
                except json.JSONDecodeError:
                    f.seek(0)
                    self.common_parameters = toml.load(f)

        if "parameters" not in self.common_parameters.keys():
            self.common_parameters = {"base-inferece": "kappa_set", "parameters": self.common_parameters}
        self.base_output = Path(output_base_path)
        analyses = []
        for lens, pars in lens_parameters.items():
            print(logger_)
            logger_.info(f"Building analysis for lens {lens}")
            analyses.append(self.build_single_analysis(lens, pars, logger, **kwargs))
        return analyses

    def build_single_analysis(self, lens_name: str, lens_paramters: dict, logger: logging.Logger, **kwargs):
        """
        Here, we just have to combine the parameters. The kappa set analysis
        will check everything to make sure its valid. Anything that is wrong with 
        the common parameters should cause this to fail on the first attempt.
        """
        system_parameters = deepcopy(self.common_parameters)
        system_parameters["parameters"].update(lens_paramters)
        system_output = self.base_output / lens_name
        system_parameters["parameters"].update({"output_base_path": system_output})
        kappa_set_template = Path(kappaSet.__file__).parents[0] / "kappa_set_template.json"
        with open(kappa_set_template, "r") as f:
            kappa_set_template_parameters = json.load(f)
        system_parameters["name"] = lens_name
        system_analysis_object = build_analysis(system_parameters, kappa_set_template_parameters, kappaSet, **kwargs)
        return system_analysis_object

class attach_wlm(Transformation):
    """
    Attaching values from the weak lensing maps to the weighted number counts
    requires user input, so we do this first to get it all done on the front end
    """
    def __call__(self, analyses: list):
        for name, analysis in analyses.items():
            print(f"Checking weak lensing maps for lens {name}")
            analysis.run_to("build_analyses")

class run_analyses(Transformation):
    def __call__(self, *args, **kwargs):
        return self.run_analyses(*args, **kwargs)
    def run_analyses(self, analyses: list):
        for name, analysis in analyses.items():
            print(f"Working on analysis for lens system {name}")
            analysis.run_analysis()

