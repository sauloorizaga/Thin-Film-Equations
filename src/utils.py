import yaml

class ConfigReader:
    def __init__(self, path) -> None:
        self.path = path
        self.data = self._read_yaml()
        self._set_simulation_attrs()
    
    def _read_yaml(self):
        with open(self.path, "r") as f:
            return yaml.safe_load(f)
        
    def _set_simulation_attrs(self):
        for key, value in self.data["simulation"].items():
            setattr(self, key, value)