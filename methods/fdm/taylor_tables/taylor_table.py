from abc import ABC, abstractmethod

class TaylorTable(ABC):

    @staticmethod
    @abstractmethod
    def lambdas():
        pass

