from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from models.model.domain import Domain
from models.model.domain import Domains
from models.model.variable import RegionEnum, Variables
from models.model.parameter import Parameters


ureg = UnitRegistry()
Q_ = ureg.Quantity


class Model(ABC):

    iter = None

    jacobian = None

    def __call__(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        return self.residue(t, y, yp, par)

    def ss_residue(self, y):

        res, _ = self.residue(0, y, np.zeros_like(y))
        print("res: {}".format(np.sum(res**2)))

        return res

    @abstractmethod
    def residue(self, t: Union[None,float], y: np.ndarray, yp: np.ndarray, par=None) -> Tuple[np.array, int]:
        pass

    class Parameters:
        pass
