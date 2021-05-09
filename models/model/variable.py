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

    def __init__(self, name, value: np.ndarray = None, unit: str = "", domains: Tuple[Domain, ...] = (), description=""):

        self.name = name
        self.domains = domains
        self.description = description
        self.shape = tuple([len(domain) for domain in domains])
        self.offset = 0

        size = 1
        for ele in  self.shape:
            size *= ele
        self.size = size

        if value is None:
            eng_value = None
            self.base_value = None
            self.base_unit = None
        else:
            eng_value = value * ureg(unit)
            self.base_value = eng_value.to_base_units().magnitude
            self.base_unit = str(eng_value.to_base_units().units)
        self.eng_value = eng_value

    def register(self, offset):
        self.offset = offset

    def parse(self, y: np.ndarray):
        ini = self.offset
        fini = ini + self.size
        return np.reshape(y[ini:fini], newshape=self.shape)

    def __call__(self):
        return self.base_value


class Variables:

    def __init__(self, variables: Tuple[Variable, ...] = ()):

        self._offset = 0
        _variables = {}
        for variable in variables:
            variable.offset = self._offset
            self._offset += variable.size
            _variables[variable.name] = variable

        self.variables = _variables

    def __get__(self) -> Dict[str, Variable]:
        return self.variables

    def __getitem__(self, key) -> Variable:
        return self.variables[key]

    def keys(self):
        return self.variables.keys()

    def read(self):
        return np.concatenate([val().flatten() for val in self.variables.values()], axis=None)

    def parse(self, y):
        return [val.parse(y) for val in self.variables.values()]