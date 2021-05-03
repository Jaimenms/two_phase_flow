from unittest import TestCase
from solvers.dassl_solver import DasslSolver
from models.burgers import Burgers
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from models.model.model_parameter import ModelParameter


class TestDasslSolverBurgers(TestCase):

    def case1(self):

        vL = 1.0
        vR = 0.5
        xi = 0.2
        N = 100

        t0 = np.linspace(0, 0.6, 2)
        xa = np.linspace(0, 1., N)
        x = np.linspace(0, 1., N)

        y0 = np.zeros(N)
        y0 = np.where(x < xi, vL, y0)
        y0 = np.where(x >= xi, vR, y0)

        lb = ModelParameter("lb", vL, "m/s")
        ub = ModelParameter("ub", vR, "m/s")

        model = Burgers(x, lb=lb, ub=ub, scheme=SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.CUBISTA2)
        model.Parameters.y_LB = vL

        return t0, y0, x, model

    def case2(self):

        N = 100

        #t0 = np.linspace(0, 0.158, 2)
        t0 = np.linspace(0, 0.5, 2)
        x = np.linspace(0, 1., N)
        y0 = +(np.sin(2*3.1415*x) + np.sin(3.1415*x)/2)

        #lb = ModelParameter("lb", 0.0, "m/s")
        #ub = ModelParameter("ub", 0.0, "m/s")
        lb = ModelParameter("lb", None, "m/s")
        ub = ModelParameter("ub", None, "m/s")

        model = Burgers(x, lb=lb, ub=ub, scheme=SchemeM1FDMEnum.CENTRAL_N4, flux_delimiter=FluxDelimiterEnum.SMART2)

        return t0, y0, x, model

    def plot_result(self, t, x, y):

        plt.figure()
        plt.plot(x, y[0,:],"k-")
        plt.plot(x, y[-1,:],"bo-")
        plt.ylabel('y')
        plt.xlabel('distance')
        plt.title('Model1 Solution')
        plt.legend(["t={}".format(t[0]), "t={}".format(t[-1])])
        plt.show()

    def test_burgers_1(self, plot=True):

        t0, y0, x, model = self.case1()

        t, y, yp = DasslSolver.run(model, t0, y0, display=True, rtol=1e-6, atol=1e-6)

        if plot:
            self.plot_result(t, x, y)
        else:
            self.assertAlmostEqual(t[-1], t0[-1])
            #self.assertAlmostEqual(y[-1][0], 0.135, places=3)

    def test_burgers_2(self, plot=True):

        t0, y0, x, model = self.case2()

        t, y, yp = DasslSolver.run(model, t0, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            self.plot_result(t, x, y)
        else:
            self.assertAlmostEqual(t[-1], t0[-1])
            #self.assertAlmostEqual(y[-1][0], 0.135, places=3)
