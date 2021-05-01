from abc import ABC, abstractmethod
import numpy as np
from scipy import sparse


class FDMWeights(ABC):

    N = 0

    lambdas = None

    @staticmethod
    @abstractmethod
    def elements(N: int, L: int, i: int):
        pass

    def __call__(self, L: int, i: int, x: np.ndarray):

        l, ini, fini = self.elements(self.N, L, i)

        x_i = x[ini:fini]
        c = np.empty(self.N+1)
        for k in range(self.N+1):
            c[k] = self.lambdas[l][k](x_i)

        return c, ini, fini

    def matrix(self, x):

        L = len(x)

        weights = np.zeros((L, L))
        for i in range(L):
            c_i, ini, fini = self.__call__(L, i, x)
            weights[i, ini:fini] = c_i

        return weights
