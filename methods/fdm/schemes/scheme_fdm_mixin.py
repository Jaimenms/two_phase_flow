from abc import ABC, abstractmethod

class SchemeFDMMixin(ABC):

    @staticmethod
    @abstractmethod
    def elements(N: int, L: int, i: int):
        pass
