from typing import List, Dict
from networkx import DiGraph, simple_cycles, isolates, draw
from lenskappa.analysis.transformation import Transformation
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import to_agraph
from networkx import draw, all_simple_paths
from pathlib import Path
import toml
from functools import reduce 
from lenskappa.locations import LENSKAPPA_CONFIG_LOCATION
import json
class InferenceException(Exception):
    pass

def build_inference(config, basic_config, module):
    """
    Build an inference from a configuration file and a module
    
    The "config" parameter should be a dictionary with the
    necessary parameters for running this inference.

    The basic_config inference should be a dictionary with
    the configuration template for this inference.

    The "module" should be the module where the transformations
    for this inference are defined.

    """
    transformations = {}
    for tname in basic_config["transformations"].keys():
        try:
            f = getattr(module, tname)
            transformations.update({tname: f()})
        except AttributeError:
            name = basic_config["name"]
            raise AttributeError(f"Unable to find transformation {tname} required by inference {name}")
    final_cfg = basic_config | config
    return Inference(transformations, final_cfg)

with open(LENSKAPPA_CONFIG_LOCATION / "base_inference.json", "r") as f:
    base_inference_config = json.load(f)
class Inference:
    """
    The inference class is a high-level class for transforming
    data. It was originally designed for computing kappa_ext
    using weighted number counts and weak lensing maps from the
    Millennium Simulation.

    The goal is to load a series of transformations, each of which
    takes in the output of other transformations and additional
    parameters as needed. "Transformation" is a very general term here.
    It could refer to simply loading data from disk, for exampl.e
    """
    base_inference_config = base_inference_config
    def __init__(self, transformations: Dict[str, Transformation], params: dict):
        self.transformations = transformations
        self.params = params
        if set(self.params["transformations"].keys()) != set(self.transformations.keys()):
                raise InferenceException("There should be one transformation in the \"transformation\""\
                                         " for each transformation in the \"transformations\" input dictonary"\
                                         " (and vice-versa).")
        self.verify_inference()
        self.internal_outputs = {}
        self.has_run = {name: False for name in self.transformations}



    def verify_inference(self):
        """
        Here, we build a dependency graph and check to make
        sure the inference is valid. For an inference to be valid, 
        it must have at least one transformation marked as an "output"
        transformation.
        
        """
        self.outputs = [k for k, v in self.params["transformations"].items() if v.get("is-output", False)]
        if not self.outputs:
            raise InferenceException("Inferences must include at least one output transformation")
        
        self.build_dependency_graph()
        self.check_params()
        self.verify_params()

    def build_dependency_graph(self):
        """
        Once an inference has been defined, we have to check its validity.
        Obvious failrue cases include if two transformations depend on the output
        of each other. Or if there is a loop (i.e. three transformations which
        depend on eachother). To check this, we consruct a directed graph of 
        transformation dependencies and check for cycles. Checking for cycles
        also ensures that a graph has at least one transformation with no dependencies,
        which is required.
        
        For each transformation, there should be an associated entry in the parameters
        dictionary which specifies dependencies. A transformation without this entry
        will be assumed to have no dependencies. This method will also halt if it finds
        an isolated transformation. That is to say, one which has no dependencies and
        depends on nothing.
        """
        transformations = self.params["transformations"]
        transformation_names = transformations.keys()
        self.dependency_graph = DiGraph()
        self.dependency_graph.add_nodes_from(transformation_names)
        self.starts = []
        for transformation, tparams in transformations.items():
            dependencies = tparams.get("dependencies", None)
            if dependencies is not None:
                if type(dependencies) != dict:
                    raise InferenceException("Dependencies should be passed as a dictionary, " \
                                            "where the key is the name of the dependency " \
                                            "transformation and the value is the name the argument " \
                                            "with its output will be assigned when passed to this "\
                                            "transformation.")
                if not all([dep in transformations.keys() for dep in dependencies]):
                    raise InferenceException("Unknown dependencies found! If this transformation needs a "\
                                             "parameter, you should put the parameter name in the needed-parameters "\
                                             "block of the dependency's configuration.")

                for dep in dependencies:
                    self.dependency_graph.add_edge(dep, transformation)
            else:
                self.starts.append(transformation)

        cycles = list(simple_cycles(self.dependency_graph))
        if cycles:
            raise InferenceException("Inference contains a dependency cycle!")
        
        isolated_transformation = list(isolates(self.dependency_graph))
        if isolated_transformation:
            raise InferenceException("Inference contains an isolated step")
        self.predecessors = {}

        # Just keep track of everything that needs to run before this transformation
        # This is used for two different functions, so we just store it. 
        for transformation in self.transformations:
            dependencies = self.dependency_graph.predecessors(transformation)
            self.predecessors.update({transformation: dependencies})

    def check_params(self):
        """
        Checks to see that all required parameters are present
        in self.parameters, and that those parameters are valid

        Throws an error if a particular parameter is not found.
        
        """
        self.needed_params = set()
        for transformation in self.transformations:
            needed_params = self.params["transformations"][transformation].get("needed-parameters", [])
            self.needed_params = self.needed_params.union(needed_params)
        self.verify_params()
        found_params = set(self.params["parameters"].keys())

        if (missing := self.needed_params - found_params):
            missing = ",".join(missing)
            raise InferenceException(f"Inference is missing required parameters {missing}")
    
    def verify_params(self):
        for param in self.needed_params:
            if param in self.base_inference_config["reserved-keys"]:
                raise InferenceException(f"Key \"{param}\" is a reserved keyword and cannot be used as"\
                                          " a parameter in an analysis definition.")



    def run_inference(self):
        """
        Runs the inference. Starts by running transformations that have
        no dependencies. Ensures that all nodes are visited exactly once.

        We've already verified that this is a valid dependency graph, i.e. that
        there is at least one transformation with no dependencies, and that there
        are no loops. We loop over every node that has not already been run, check
        to see if its dependencies HAVE been run, and run it if so. Once everything
        has been run once, we break.
        """

        for transformation in self.starts:
            self.run_transformation(transformation)
        while True:
            for transformation in self.transformations:
                if self.has_run[transformation]:
                    continue
                dependencies = self.predecessors[transformation]
                if all([self.has_run[k] for k in dependencies]):
                    self.run_transformation(transformation)
            if all(self.has_run.values()):
                break
        return self.cleanup()


    def run_transformation(self, name):
        """
        Runs a single transformation by name. Before we get here,
        we have always checked to make sure all the necessary parameters are
        present. We also make sure to run all the transformations in the correct
        order, so any transformations this transformation depends on should have
        been run and their output should have been cached. As such, this function
        should pretty much never fail, unless the underlying transformation
        fails. 
        
        """
        print(f"MAIN: Running transformation \"{name}\"")

        arguments = {}
        for dep in self.dependency_graph.predecessors(name):
            alias = self.params["transformations"][name]["dependencies"][dep]
            if not alias:
                continue
            output = self.internal_outputs[dep]
            arguments.update({alias: output})
        needed_params = self.params["transformations"][name].get("needed-parameters", [])
        optional_params = self.params["transformations"][name].get("optional-parameters", [])
        all_params = needed_params + optional_params
        for param in all_params:
            pvalue = self.params["parameters"].get(param, False)
            if pvalue:
                arguments.update({param: pvalue})
        transformation_output = self.transformations[name](**arguments)
        self.internal_outputs.update({name: transformation_output})
        self.has_run[name] = True

    def run_to(self, transformation_name: str):
        """
        Run this analysis up to and including a particular transformation. This will
        determine the shortest path that includes all prerequisite transformations.
        In other words, it will do the least amount of work possible to get to this 
        transformation. If you then call "run," the analysis will pick up
        where it left off. Alternatively, call "run_to" again to advance the analysis
        in steps.
        """
        if transformation_name not in self.transformations:
            raise KeyError(f"Transformation {transformation_name} not found in this analysis!")

        #If this transformation has no prereqs, just run it
        if transformation_name in self.starts:
            self.run_transformation(transformation_name)
            return self.internal_outputs[transformation_name]
        
        #If this transformation is the only output transformation, throw an error
        if transformation_name in self.outputs and len(self.outputs) == 0:
            raise InferenceException(f"Running to transformation {transformation_name} would require running"\
                                    "the full analysis! Use \"run_analysis\" instead.")
        
        #Otherwise, we determine which transformations we need to run first,
        #and run those
        necessary_paths = []
        for start in self.starts:
            necessary_paths.extend(list(all_simple_paths(self.dependency_graph, start, transformation_name)))
        #Get the nodes in these paths that are actually unique
        unique_transformations = reduce(lambda l, r: set(l).union(r), necessary_paths )
        while True:
            for transformation in unique_transformations:
                if self.has_run[transformation]:
                    continue
                dependencies = self.predecessors[transformation]
                if all([self.has_run[k] for k in dependencies]):
                    self.run_transformation(transformation)
            if self.has_run[transformation_name]:
                break
        return self.internal_outputs[transformation_name]

    def reset(self):
        """
        Completely reset an inference, as though nothing has run. Does
        not discard the dependency graph and related products.
        This cannot be undone, so use carefully.

        My use case involved a transformation that pulls data from a secondary
        source and saves it to one of the other input files, but requires user
        input to do so. If I would like to run many of these at once, I want
        the user to be able to provide these inputs all at the beginning of the
        analysis, rather than periodically throughout the (potentially very long)
        runtime. The transformation checks to see if this copy has already happened,
        so we reset and run the whole analysis a second time when no input is required.
        """
        self.internal_outputs = {}
        self.has_run = {name: False for name in self.transformations}


    def cleanup(self):
        if len(self.outputs) == 1:
            return self.internal_outputs[self.outputs[0]]
        else:
            return {name: self.internal_outputs[name] for name in self.outputs}

    
