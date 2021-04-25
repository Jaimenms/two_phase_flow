from abc import ABC, abstractmethod


class FluxDelimiter(ABC):

    @abstractmethod
    def __call__(self, phi_p: float, tetha_f: float, tetha_p):
        return phi_p