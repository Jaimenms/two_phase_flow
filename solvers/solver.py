from abc import ABC, abstractmethod
from models.model.model import Model
import numpy as np

class Solver(ABC):

    @staticmethod
    @abstractmethod
    def run(self, model: Model, t: np.ndarray, y0: np.ndarray, yp0=None, rtol=1e-6, atol=1e-8, index=None):
        pass
