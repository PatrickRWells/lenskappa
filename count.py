import parser
import numpy as np
import weightfns

class count:
    def __init__(self, name, config):
        self.name = name
        self.eqn = self._parse_equation()

    def _parse_conditions(self, param):
        return True
    
    def _parse_equation(self):
        return getattr(weightfns, self.name)


if __name__ == '__main__':
    import toml
    config = toml.load('config/counts.toml')
    weight = count('zweight', config)