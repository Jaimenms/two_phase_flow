import numpy as np


class Dimensionless:

    @staticmethod
    def reynolds( D, v, rho, mu):
        """

        :param D: viameter
        :param v: velocity
        :param rho: density
        :param mu: viscosity
        :return:
        """
        return D * np.abs(v) * rho / mu