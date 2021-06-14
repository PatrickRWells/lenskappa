import toml

class ConfigManager:
    def __init__(self, config):
        self.top_config = toml.load(config)
        self.setup()
    def setup(self):
        pass

    def get_submodule_config(self, object):
        pass