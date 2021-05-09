from scipy.optimize import root_scalar
import numpy as np

class Geometry:

    @staticmethod
    def area(D):
        return 0.25 * np.pi * D**2

    @staticmethod
    def stratified_angle_res(betha, alphaL):
        return betha - np.sin(betha) - 2 * np.pi * alphaL

    @staticmethod
    def stratified_angle_fprime(betha, *args):
        return 1 - np.cos(betha)

    @staticmethod
    def stratified_angle_fprime2(betha, *args):
        return np.sin(betha)

    @staticmethod
    def stratified_angle_approximate(alphaL):
        """
        Keplers solution
        :param alphaL: liquid area fraction
        :return:
        """

        # Keplers solution
        n = 10
        M = 2 * np.pi * alphaL
        betha = M
        for i in range(n):
            betha = M + np.sin(betha)

        return betha


    @staticmethod
    def stratified_angle(alphaL):

        betha = Geometry.stratified_angle_approximate(alphaL)

        if not isinstance(alphaL, np.ndarray):

            sol = root_scalar(
                Geometry.stratified_angle_res,
                x0=betha,
                fprime=Geometry.stratified_angle_fprime,
                fprime2=Geometry.stratified_angle_fprime2,
                bracket=[0, 2*np.pi],
                method='newton',
                args=(alphaL,)
            )
            betha = sol.root

        else:
            N = len(betha)
            for i in range(N):
                sol = root_scalar(
                    Geometry.stratified_angle_res,
                    x0=betha[i],
                    fprime=Geometry.stratified_angle_fprime,
                    fprime2=Geometry.stratified_angle_fprime2,
                    bracket=[0, 2 * np.pi],
                    method='newton',
                    args=(alphaL[i],)
                )
                betha[i] = sol.root

        return betha


    @staticmethod
    def stratified_perimeters(D, betha):
        """

        :param D: diamater
        :param betha: stratified angle
        :return: tuple with liquid, gas and interface perimiters
        """

        perL = 0.5 * betha * D
        perG = np.pi * D - perL
        perI = D * np.sin(0.5 * betha)

        return perL, perG, perI
