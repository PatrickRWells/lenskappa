from lenskappa.datasets.dataset import SkyDataSet

class Simulation(SkyDataSet):
    
    def __init___(self, name, *args, **kwargs):
        """
        A simulation differs from a survey in that it doesn't
        actually represent a physical area on the sky. However, it should be able to
        be projected onto a fake sky.
        """
        super().__init__(name, *args, **kwargs)