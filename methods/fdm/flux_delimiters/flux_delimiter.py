from abc import ABC, abstractmethod
import numpy as np


class FluxDelimiter(ABC):

    @abstractmethod
    def __call__(self, phi_p: np.ndarray, tetha_f: np.ndarray, tetha_p: np.ndarray) -> np.ndarray:
        return phi_p