from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from models.model.model_domain import ModelDomain
from models.model.model_variable import ModelVariable, ModelRegionEnum
from models.model.model_parameter import ModelParameter
from models.model.model_constant import ModelConstant

ureg = UnitRegistry()
Q_ = ureg.Quantity


class Model(ABC):

    iter = None

    def __call__(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        if self.iter is None:
            self.iter = 0
        else:
            new_iter = np.count_nonzero(par.solver_t <= t)
            if new_iter > self.iter:
                self.iter = new_iter
                print("{}/{}".format(self.iter, len(par.solver_t)))

        return self.residue(t, y, yp, par)

    @abstractmethod
    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None) -> Tuple[np.array, int]:
        pass

    def jacobian(self, t: float, y: np.ndarray, yp: np.ndarray, par=None) -> Union[None, Tuple[np.array, int]]:
        return None

    @staticmethod
    @abstractmethod
    def str_equation(self) -> str:
        pass

    class Parameters:
        pass


    def __init__(self):
        self.domains = dict()
        self.constants = dict()
        self.parameters = dict()
        self.variables = dict()
        self._offset = 0

    def register_domain(self, input: ModelDomain):
        self.domains[input.name] = input

    def register_domains(self, input: Tuple[ModelDomain,...]):
        for val in input:
            self.register_domain(val)

    def register_variable(self, input: ModelVariable):
        input.offset = self._offset
        self._offset += input.size
        self.variables[input.name] = input

    def register_variables(self, input: Tuple[ModelVariable,...]):
        for val in input:
            self.register_variable(val)

    def register_parameter(self, input: ModelParameter):
        self.parameters[input.name] = input

    def register_parameters(self, input: Tuple[ModelParameter,...]):
        for val in input:
            self.register_parameter(val)

    def register_constant(self, input: ModelConstant):
        self.constants[input.name] = input

    def register_constants(self, input: Tuple[ModelConstant,...]):
        for val in input:
            self.register_constant(val)

    def parse(self, name, y: np.ndarray):
        variable = self.variables[name]
        ini = variable.offset
        fini = variable.offset + variable.size
        return np.reshape(y[ini:fini], newshape=variable.shape)

    @staticmethod
    def apply_regions(value: np.ndarray, regions: Tuple[ModelRegionEnum, ...]):
        slices = tuple([region.value for region in regions])
        return value[slices]
