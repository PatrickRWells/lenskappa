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
import logging
import dask
import uuid
from dask.distributed import Client, get_client
from loguru import logger
import sys
class AnalysisException(Exception):
    pass

def build_analysis(config, basic_config, module, **kwargs):
    """
    Build an analysis from a configuration file and a module
    
    The "config" parameter should be a dictionary with the
    necessary parameters for running this analysis.

    The basic_config analysis should be a dictionary with
    the configuration template for this analysis.

    The "module" should be the module where the transformations
    for this analysis are defined.

    """
    transformations = {}
    for tname in basic_config["transformations"].keys():
        try:
            f = getattr(module, tname)
            transformations.update({tname: f()})
        except AttributeError:
            name = basic_config["name"]
            raise AttributeError(f"Unable to find transformation {tname} required by analysis {name}")
    final_cfg = basic_config | config
    
    if (setup := getattr(module, "setup", False)):
        extra_parameters = setup(final_cfg)
        final_cfg["parameters"].update(extra_parameters)

    return Analysis(transformations, final_cfg, **kwargs)

with open(LENSKAPPA_CONFIG_LOCATION / "base_analysis.json", "r") as f:
    base_analysis_config = json.load(f)
class Analysis:
    """
    The analysis class is a high-level class for transforming
    data. It was originally designed for computing kappa_ext
    using weighted number counts and weak lensing maps from the
    Millennium Simulation.

    The goal is to load a series of transformations, each of which
    takes in the output of other transformations and additional
    parameters as needed. "Transformation" is a very general term here.
    It could refer to simply loading data from disk, for example.
    """
    base_analysis_config = base_analysis_config
    def __init__(self, transformations: Dict[str, Transformation], params: dict, **kwargs):
        self.transformations = transformations
        self.params = params
        self.tags = self.params.get("tags", [])
        self.name = self.params["name"]

        self.get_clients(**kwargs)

        if set(self.params["transformations"].keys()) != set(self.transformations.keys()):
                raise AnalysisException("There should be one transformation in the \"transformation\""\
                                         " for each transformation in the \"transformations\" input dictonary"\
                                         " (and vice-versa).")

        self.verify_analysis()

        self.build_dask_graph()

        self.internal_outputs = {}
        self.has_run = {name: False for name in self.transformations}

    def get_clients(self, sub_analysis = False, logger_ = None, **kwargs):
        self.logger = logger_
        if self.logger is None:
            self.logger = logger
            self.logger.remove()
            self.logger.add(sys.stderr, level="INFO")

        if sub_analysis:
            self.client = get_client()
        else:
            nthreads = self.params.get("threads", 1)
            self.logger.info(f"Found threads = {nthreads}. Spinning up worker processes...")
            self.client = Client(n_workers = nthreads)

    def verify_analysis(self):
        """
        Here, we build a dependency graph and check to make
        sure the analysis is valid. For an analysis to be valid, 
        it must have at least one transformation marked as an "output"
        transformation.
        
        """
        self.outputs = [k for k, v in self.params["transformations"].items() if v.get("is-output", False)]
        if not (self.outputs or "analysis-builder" in self.tags) :
            raise AnalysisException("Analysis must include at least one output transformation")
        
        self.build_dependency_graph()
        self.check_params()
        self.verify_params()

    def build_dependency_graph(self):
        """
        Once an analysis has been defined, we have to check its validity.
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
                    raise AnalysisException("Dependencies should be passed as a dictionary, " \
                                            "where the key is the name of the dependency " \
                                            "transformation and the value is the name the argument " \
                                            "with its output will be assigned when passed to this "\
                                            "transformation.")
                if not all([dep in transformations.keys() for dep in dependencies]):
                    raise AnalysisException("Unknown dependencies found! If this transformation needs a "\
                                             "parameter, you should put the parameter name in the needed-parameters "\
                                             "block of the dependency's configuration.")

                for dep in dependencies:
                    self.dependency_graph.add_edge(dep, transformation)
            else:
                self.starts.append(transformation)

        cycles = list(simple_cycles(self.dependency_graph))
        if cycles:
            raise AnalysisException(f"Analysis {self.name} contains a dependency cycle!")
        
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
            raise AnalysisException(f"Analysis is missing required parameters {missing}")
    
    def build_dask_graph(self, *args, **kwargs):
        self.delayed = {}
        self.scheduled = {n: False for n in self.transformations}
        while True:
            for transformation in self.transformations:
                if self.scheduled[transformation]:
                    continue
                dependencies = self.predecessors[transformation]
                if not dependencies or all([self.scheduled[k] for k in dependencies]):
                    self.schedule_transformation(transformation)
            if all(self.scheduled.values()):
                break

    def schedule_transformation(self, transformation):
        if self.params["transformations"][transformation].get("analysis-builder", False):
            self.handle_analysis_builder(transformation)
        else:
            inputs = self.get_transformation_parameters(transformation, scheduled = True)
            delayed_fn = dask.delayed(self.transformations[transformation])(**inputs)
            self.delayed.update({transformation: delayed_fn})
            self.scheduled[transformation] = True
    
    def handle_analysis_builder(self, transformation):
        """
        When one of the transformations in an analysis creates new analyses,
        we have to be careful when adding them to dask. Ideally, we'd be able
        to fully expand the graph before running anything. But if this analysis
        builder depends on results of pervious transformations, there's no way
        to avoid running this in bits and pieces.  
        
        So we run any transformation leading up to this builder, then take the 
        delatedy outputs from the analyses that have been built and glue it onto 
        the execution graph we've already built. 
        """
        
        predecessors = self.predecessors[transformation]
        for predecessor in predecessors:
            self.run_to(predecessor)
        
        analyses = self.run_single(transformation, sub_analysis = True)
        self_is_output = self.params["transformations"][transformation].get("is-output", False)
        for a in analyses:
            outputs = [a.delayed[output] for output in a.outputs]
            if self_is_output:
                a_uuid = str(uuid.uuid1())
                output_names = ["-".join([output, a_uuid]) for output in a.outputs]
                output_dict = dict(list(zip(output_names, outputs)))
                self.delayed.update(output_dict)
                self.outputs.extend(output_names)
            else:
                self.delayed.update({transformation: outputs})
                
        self.outputs.remove(transformation)
        self.scheduled[transformation] = True

    

    def verify_params(self):
        for param in self.needed_params:
            if param in self.base_analysis_config["reserved-keys"]:
                raise AnalysisException(f"Key \"{param}\" is a reserved keyword and cannot be used as"\
                                          " a parameter in an analysis definition.")



    def run_analysis(self):
        """
        Runs the analysis. Starts by running transformations that have
        no dependencies. Ensures that all nodes are visited exactly once.

        We've already verified that this is a valid dependency graph, i.e. that
        there is at least one transformation with no dependencies, and that there
        are no loops. We loop over every node that has not already been run, check
        to see if its dependencies HAVE been run, and run it if so. Once everything
        has been run once, we break.
        """
        outputs = {}
        for output in self.outputs:
            outputs.update({output: self.delayed[output].compute()})
        return outputs



    def get_transformation_parameters(self, name, scheduled = False):
        """
        Retrieve any parameters that are required by the transformation with this name
        and return them as a dictionary. 
        
        """

        arguments = {}
        for dep in self.dependency_graph.predecessors(name):
            alias = self.params["transformations"][name]["dependencies"][dep]
            if not alias:
                continue
            if scheduled:
                output = self.delayed[dep]
            else:
                output = self.internal_outputs[dep]
            arguments.update({alias: output})
        needed_params = self.params["transformations"][name].get("needed-parameters", [])
        optional_params = self.params["transformations"][name].get("optional-parameters", [])
        all_params = needed_params + optional_params
        for param in all_params:
            pvalue = self.params["parameters"].get(param, None)
            if pvalue is not None: #ugly I know, but "if val" doesn't work for some values
                arguments.update({param: pvalue})
        arguments.update({"logger_": self.logger})
        arguments.update({"name": self.name})
        return arguments

    def run_single(self, name, **kwargs):
        """
        Runs a single transformation by name. Before we get here,
        we have always checked to make sure all the necessary parameters are
        present. We also make sure to run all the transformations in the correct
        order, so any transformations this transformation depends on should have
        been run and their output should have been cached. As such, this function
        should pretty much never fail, unless the underlying transformation
        fails. 
        
        """
        self.logger.debug(f"Running transformation \"{name}\"")
        transformation_parameters = self.get_transformation_parameters(name, )
        transformation_output = self.transformations[name](**transformation_parameters, **kwargs)
        return transformation_output

    def run_to(self, transformation_name: str, **kwargs):
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
        
        #If this transformation is the only output transformation, throw an error
        
        self.run_transformation(transformation_name, **kwargs)
        #Otherwise, we determine which transformations we need to run first,
        #and run those

    def reset(self):
        """
        Completely reset an analysis, as though nothing has run. Does
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

    
