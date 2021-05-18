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

    @staticmethod
    def floud(rhoL, rhoG, D, vLS, vGS):
        """

        :param rhoL:
        :param rhoG:
        :param D:
        :param vLS:
        :param vGS:
        :return:
        """

        FL = (rhoL / (rhoL - rhoG) / 9.81 / D) ** 0.5 * np.abs(vLS)
        FG = (rhoG / (rhoL - rhoG) / 9.81 / D) ** 0.5 * np.abs(vGS)

        return FL, FG