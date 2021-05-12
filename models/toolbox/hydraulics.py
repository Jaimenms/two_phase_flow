import numpy as np

class Hydraulics:

    @staticmethod
    def ff_via_churchill(Re, epw, D):
        """

        :param Re: reynolds number
        :param Epw: pipe roughness
        :param D: Diameter
        :return: tuple with friction factor and Darcy factor
        """

        A = (2.457 * np.log(((7 / Re) ** 0.9 + 0.27 * epw / D) ** -1)) ** 16
        B = (37530 / Re) ** 16
        ff = 2 * ((8 / Re) ** 12 + (A + B) ** -1.5) ** (1 / 12)
        fD = 4 * ff

        return ff, fD

    @staticmethod
    def shear_stress(fD, rho, v):
        """

        :param fD: Darcy factor
        :param rho: density
        :param v: velocity
        :return:
        """

        return 0.5 * fD * rho * v * np.abs(v)

    @staticmethod
    def velocity(q, rho, A):
        """

        :param q: mass flowrate
        :param rho: density
        :param A: area
        :return:
        """
        return q / rho / A