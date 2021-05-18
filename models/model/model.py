from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry


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

    def numerical_jacobian_y_block(self, t, y, yp, par):
        n = len(y)
        jacAs = []
        for i in range(n):
            y1 = 1*y
            y1[i] = y[i] + 1e-6*np.abs(y[i]) + 1e-9
            res1, _ = self.residue(t=t, y=y1, yp=yp, par=par)
            y2 = 1*y
            y2[i] = y[i] - 1e-6*np.abs(y[i]) - 1e-9
            res2, _ = self.residue(t=t, y=y2, yp=yp, par=par)
            jacA_i = (res1 - res2)/(y1[i] - y2[i])
            jacAs.append(jacA_i)
        return np.vstack(jacAs)

    def numerical_jacobian_y_block_vectorized(self, t, y, yp, par):
        n = len(y)

        Y1 = []
        Y2 = []
        YP = []

        DEN = []
        for i in range(n):
            y1 = 1*y
            y1[i] = y[i] + 1e-6*np.abs(y[i]) + 1e-9
            Y1.append(y1)
            y2 = 1*y
            y2[i] = y[i] - 1e-6*np.abs(y[i]) - 1e-9
            Y2.append(y2)
            YP.append(yp)
            DEN.append((y1[i]-y2[i])*np.ones_like(y1))

        Y1 = np.vstack(Y1).T
        Y2 = np.vstack(Y2).T
        YP = np.vstack(YP).T
        DEN  = np.vstack(DEN).T

        res1, _ = self.residue(t=t, y=Y1, yp=YP, par=par)
        res2, _ = self.residue(t=t, y=Y2, yp=YP, par=par)

        NUM =np.reshape(res1 - res2,(n,n))
        out = NUM/DEN
        return out.T

    def numerical_jacobian_yprime_block(self, t, y, yp, par):
        n = len(y)
        jacBs = []
        for i in range(n):
            y1p = 1*yp
            y1p[i] = yp[i] + 1e-6*np.abs(yp[i]) + 1e-9
            res1p, _ = self.residue(t=t, y=y, yp=y1p, par=par)
            y2p = 1*yp
            y2p[i] = yp[i] - 1e-6*np.abs(yp[i]) - 1e-9
            res2p, _ = self.residue(t=t, y=y, yp=y2p, par=par)
            jacB_i = (res1p - res2p)/(y1p[i] - y2p[i])
            jacBs.append(jacB_i)
        out = np.vstack(jacBs)
        return out

    def numerical_jacobian(self, t, y, yp, cj, par):

        print("Jac calc at {} s".format(t))

        jacA = self.numerical_jacobian_y_block_vectorized(t, y, yp, par)
        jacB = self.numerical_jacobian_yprime_block(t, y, yp, par)

        jac = jacA + cj * jacB

        return jac.T, 0

    def numerical_jacobian_2(self, t, y, yp, cj, par):

        print("Jac calc at {} s".format(t))

        jacA = self.numerical_jacobian_y_block(t, y, yp, par)
        jacB = self.numerical_jacobian_yprime_block(t, y, yp, par)

        jac = jacA + cj * jacB

        return jac.T, 0

    class Parameters:
        pass
