from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from models.model.domain import Domain
from models.model.variable import RegionEnum
import matplotlib.pyplot as plt

ureg = UnitRegistry()
Q_ = ureg.Quantity


class Model(ABC):

    iter = None

    def __call__(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        return self.residue(t, y, yp, par)

    def ss_residue(self, y):

        res, _ = self.residue(0, y, np.zeros_like(y))
        print("res: {}".format(np.sum(res**2)))

        return res

    @abstractmethod
    def residue(self, t: Union[None,float], y: np.ndarray, yp: np.ndarray, par=None) -> Tuple[np.array, int]:
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
        self.parameters = dict()
        self.variables = dict()

    def register_domain(self, input: Domain):
        self.domains[input.name] = input

    def register_domains(self, input: Tuple[Domain, ...]):
        for val in input:
            self.register_domain(val)

    @staticmethod
    def apply_regions(value: np.ndarray, regions: Tuple[RegionEnum, ...]):
        slices = tuple([region.value for region in regions])
        return value[slices]

    def plot_result(self, t, y):

        x = self.domains["x"]()

        y0 = y[0]

        for i, ti in enumerate(t):
            if i > 0:
                yi = y[i]
                keys = list(self.variables.keys())
                fig, axs = plt.subplots(len(keys))
                for j, var in enumerate(keys):
                    yj = self.variables[var].parse(yi)
                    y0j = self.variables[var].parse(y0)
                    axs[j].plot(x, y0j, "k-")
                    axs[j].plot(x, yj, "bo-")
                    axs[j].set_ylabel(var)
                    axs[j].set_xlabel('distance')
                    axs[j].legend(["t={}".format(t[0]), "t={}".format(t[i])])

        plt.show()
