import numpy as np
from typing import Tuple, Union, Dict
from pint import UnitRegistry
from models.model.domain import Domain, Domains
from enum import Enum, unique

ureg = UnitRegistry()
Q_ = ureg.Quantity


@unique
class RegionEnum(Enum):
    CLOSED_CLOSED = np.s_[0:]
    CLOSED_OPEN = np.s_[0:-1]
    OPEN_OPEN = np.s_[1:-1]
    OPEN_CLOSED = np.s_[1:]
    LOWER = np.s_[0]
    UPPER = np.s_[-1]
    ALL = np.s_[:]
    LOWER_PLUS_ONE = np.s_[1]
    UPPER_MINUS_ONE = np.s_[-2]



class Variable:

    def __init__(self, name, unit: str = "", domains: Tuple[Domain, ...] = (), description=""):

        self.name = name
        self.domains = domains
        self.description = description
        self._value = None
        self._unit = str(ureg(unit).to_base_units().units)
        self.shape = tuple([len(domain) for domain in domains])
        self.offset = 0
        self.indices = []

        size = 1
        for ele in  self.shape:
            size *= ele
        self.size = size

    def register(self, offset):
        self.offset = offset
        self.indices = np.array(list(range(offset, offset+self.size)), dtype="int64")


    @property
    def value(self):
        if self._value is None:
            raise ValueError("Sorry but the parameter value is not defined.")
        return self._value

    def set_ic(self, value: Union[float, np.ndarray], unit: str= ""):

        value_base_unit = str(ureg(unit).to_base_units().units)

        if self._unit != value_base_unit:
            raise ValueError("Sorry but the parameter unit is not consistent")
        self._value = value * ureg(unit).to_base_units().magnitude

    def parse(self, y: np.ndarray) -> np.ndarray:
        yi = np.take(y, self.indices, axis=0)
        if yi.ndim == 1:
            shape = self.shape
        else:
            shape = (*self.shape, *yi.shape[1:])
        return np.reshape(yi, newshape=shape)


    def __call__(self) -> np.ndarray:
        return self._value


class Variables:

    def __init__(self, variables: Tuple[Variable, ...] = ()):

        self._offset = 0
        _variables = {}
        for variable in variables:
            variable.register(self._offset)
            self._offset += variable.size
            _variables[variable.name] = variable

        self.variables = _variables

    def __getitem__(self, key) -> Variable:
        return self.variables[key]

    def keys(self):
        return self.variables.keys()

    def get_ic_array(self):
        return np.concatenate([val().flatten() for val in self.variables.values()], axis=None)

    def parse(self, y):
        return [val.parse(y) for val in self.variables.values()]