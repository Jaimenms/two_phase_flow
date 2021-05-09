from typing import Union, Tuple, Dict
from pint import UnitRegistry
from scipy import interpolate
import numpy as np
from abc import ABC, abstractmethod

ureg = UnitRegistry()


class Parameter(ABC):

    def __init__(self, name: str, value: Union[float, np.ndarray]=None, unit: str="", tspan= None, description=""):

        if value is None:
            eng_value = None
            self.base_value = None
            self.base_unit = None
        else:
            eng_value = value * ureg(unit)
            self.base_value = eng_value.to_base_units().magnitude
            self.base_unit = str(eng_value.to_base_units().units)
        self.name = name
        self.eng_value = eng_value
        self.description = description

    @abstractmethod
    def __call__(self):
        pass


class TimeDependentParameter(Parameter):

    def __init__(self, name: str, value: Union[float, np.ndarray], unit: str, description="", tspan= None):
        super().__init__(name, value, unit, description, tspan)
        self.f = interpolate.interp1d(tspan, self.base_value)

    def __call__(self, t):
        return self.f(t)


class ConstantParameter(Parameter):

    def __init__(self, name: str, value: Union[float, np.ndarray], unit: str, description=""):
        super().__init__(name, value, unit, description)

    def __call__(self):
        return self.base_value


class Parameters:

    def __init__(self, parameters: Tuple[Parameter, ...] = ()):
        self._offset = 0
        _parameters = {}
        for parameter in parameters:
            _parameters[parameter.name] = parameter

        self.parameters = _parameters

    def __get__(self) -> Dict[str, Parameter]:
        return self.parameters

    def __getitem__(self, key) -> Parameter:
        return self.parameters[key]


