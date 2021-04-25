from methods.fdm.fdm_mixin import FDMMixin
from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple


class Model(ABC, FDMMixin):

    def __init__(self):
        pass

    @abstractmethod
    def residue(self, t: float, y: np.ndarray, yp: np.ndarray) -> Tuple[np.array, int]:
        pass

    @staticmethod
    @abstractmethod
    def str_equation(self) -> str:
        pass
