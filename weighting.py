# P R Wells
# This code is based on earlier code by CE Rusu. It is used to calculate weights for objects in a lens field. 

from lens import Lens
import toml


def init_lens(name):
    return Lens(name)
def load_constants():
    config_location = "config/constants.toml"
    return toml.load(config_location)


if __name__ == "__main__":
    init_lens("HE1104")
    data = load_constants()
    print(data)
