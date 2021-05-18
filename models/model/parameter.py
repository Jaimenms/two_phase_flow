from typing import Union, Tuple, Dict
from pint import UnitRegistry
from scipy import interpolate
import numpy as np
from abc import ABC, abstractmethod

ureg = UnitRegistry()


class Parameter(ABC):

    def __init__(self, name: str, unit: str="", description=""):
        self._value = None
        self._unit = str(ureg(unit).to_base_units().units)
        self.name = name
        self.description = description

    @abstractmethod
    def __call__(self, t=None):
        pass


class TimeDependentParameter(Parameter):

    def __init__(self, name: str, unit: str, description=""):
        super().__init__(name, unit, description)

    @property
    def value(self):
        if self._value is None:
            raise ValueError("Sorry but the parameter value is not defined.")
        return self._value

    def set(self, time_span: np.ndarray, value_span: np.ndarray, time_unit: str="", value_unit: str=""):

        value_base_unit = str(ureg(value_unit).to_base_units().units)

        if self._unit != value_base_unit:
            raise ValueError("Sorry but the parameter unit is not consistent")

        time_span = time_span * ureg(time_unit).to_base_units().magnitude
        value_span = value_span * ureg(value_unit).to_base_units().magnitude

        self._value = interpolate.interp1d(time_span, value_span)

    def __call__(self, t=None) -> Union[float, np.ndarray]:
        return self._value(t)


class ConstantParameter(Parameter):

    def __init__(self, name: str, unit: str, description=""):
        super().__init__(name, unit, description)

    @property
    def value(self):
        if self._value is None:
            raise ValueError("Sorry but the parameter value is not defined.")
        return self._value

    def set(self, new_value: Union[float, np.ndarray], value_unit: str= ""):

        value_base_unit = str(ureg(value_unit).to_base_units().units)

        if self._unit != value_base_unit:
            raise ValueError("Sorry but the parameter unit is not consistent")

        if new_value is None:
            self._value = None
        else:
            self._value = new_value * ureg(value_unit).to_base_units().magnitude

    def __get__(self):
        return self._value

    def __call__(self, t=None) -> Union[float, np.ndarray]:
        return self._value


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

    def keys(self):
        return self.parameters.keys()


