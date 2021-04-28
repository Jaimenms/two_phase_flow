import numpy as np
from pint import UnitRegistry

ureg = UnitRegistry()


class ModelDomain:

    def __init__(self, name, value: np.ndarray = None, unit="", description=""):

        if not ureg(unit).check('[length]'):
            raise TypeError("Domain unit does not have length dimension.")

        eng_value = value * ureg(unit)
        self.name = name
        self.eng_value = eng_value
        self.base_value = eng_value.to_base_units().magnitude
        self.base_unit = str(eng_value.to_base_units().units)
        self.description = description

    def __get__(self):
        return self.base_value
