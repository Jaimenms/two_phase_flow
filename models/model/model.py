from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry

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
