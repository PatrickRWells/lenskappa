from abc import ABC, abstractmethod, abstractproperty
from typing import final



class LenskappaAnalysis(ABC):

    def __init__(self, *args, **kwargs):
        """
        The LenskappaAnalysis is a core class which defines
        how a particular analysis will be performed. Every analysis
        needs to have defined inputs, and a defined output handler (or several,
        depending on the use case). It should also define a "get_step" method,
        which takes the inputs, performs all the logic that needs to be performed,
        and then returns an object that can be interpreted by the output
        handler.

        Any configuration needed by the analysis should be passed to __init__.
        For more complex anlyses with many parameters, it may be worthwhile
        to define factory functions to make this step easier.
        """

    @abstractmethod
    def get_step(self, *args, **kwargs):
        """
        Performs a single step in the analysis.        
        """
        pass
    
    @final
    def verify_setup(self, *args, **kwargs):
        """
        Automatically called when the analysis is run
        """
        pass

    @final
    def run_analysis(self, *args, **kwargs):
        """
        Runs the analysis.
        """
        pass
    