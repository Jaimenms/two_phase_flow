from unittest import TestCase
from solvers.dassl_solver import DasslSolver
from models.burgers import Burgers
import numpy as np
import matplotlib.pyplot as plt
from methods.fdm.fdm_mixin import FluxDelimiterEnum, FDMEnum

class TestDasslSolverBurgers(TestCase):

    def case1(self):

        vL = 1.0
        vR = 0.5
        xi = 0.2
        N = 200

        t0 = np.linspace(0, 0.5, 2)
        x = np.linspace(0, 1., N)

        y0 = np.zeros(N)
        y0 = np.where(x < xi, vL, y0)
        y0 = np.where(x >= xi, vR, y0)

        model = Burgers(x, order=FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.MINMOD)

        return t0, y0, x, model

    def case2(self):

        N = 150

        t0 = np.linspace(0, 0.2, 2)
        x = np.linspace(0, 1., N)
        y0 = np.sin(2*3.1415*x) + np.sin(3.1415*x)/2
        model = Burgers(x, order=FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.MINMOD)

        return t0, y0, x, model

    def plot_result(self, t, x, y):

        plt.figure()
        plt.plot(x, y[0,:],"k-")
        plt.plot(x, y[-1,:],"b-")
        plt.ylabel('y')
        plt.xlabel('distance')
        plt.title('Model1 Solution')
        plt.legend(["t={}".format(t[0]), "t={}".format(t[-1])])
        plt.show()

    def test_burgers_1(self, plot=True):

        t0, y0, x, model = self.case1()

        t, y, yp = DasslSolver.run(model, t0, y0, None, 1e-6, 1e-8)

        if plot:
            self.plot_result(t, x, y)
        else:
            self.assertAlmostEqual(t[-1], t0[-1])
            #self.assertAlmostEqual(y[-1][0], 0.135, places=3)

    def test_burgers_2(self, plot=True):

        t0, y0, x, model = self.case2()

        t, y, yp = DasslSolver.run(model, t0, y0, None, 1e-6, 1e-8)

        if plot:
            self.plot_result(t, x, y)
        else:
            self.assertAlmostEqual(t[-1], t0[-1])
            #self.assertAlmostEqual(y[-1][0], 0.135, places=3)