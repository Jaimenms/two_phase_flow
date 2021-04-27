from methods.fdm.fdm_mixin import FDMMixin
from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union

class Model(ABC, FDMMixin):

    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, t: float, y: np.ndarray, yp: np.ndarray, par=None) -> Tuple[np.array, int]:
        pass

    def jacobian(self, t: float, y: np.ndarray, yp: np.ndarray, par=None) -> Union[None, Tuple[np.array, int]]:
        return None

    @staticmethod
    @abstractmethod
    def str_equation(self) -> str:
        pass

    class Parameters:
        pass