import numpy as np
from itertools import product
from abc import ABC, abstractmethod


class Operation(ABC):

    axis = 0

    @abstractmethod
    def slice_call(self, f: np.ndarray, *args_slice):
        pass

    def recursive_slice_call(self, f, *args):

        Ni, Nk = f.shape[:self.axis], f.shape[self.axis + 1:]

        grads = np.empty_like(f)
        for elem in product(np.ndindex(Ni), np.ndindex(Nk)):
            ii, kk = elem
            ind = ii + np.s_[:, ] + kk
            f_slice = f[ind]
            args_slice = tuple(arg[ind] for arg in args)
            grads[ind] = self.slice_call(f_slice, *args_slice)

        return grads