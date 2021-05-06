from unittest import TestCase
from methods.fdm.adaptative_grid import AdaptativeGrid
import numpy as np
from models.model.model_parameter import ModelParameter
from solvers.dassl_solver import DasslSolver
from models.burgers import Burgers
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
import matplotlib.pyplot as plt
import time


class TestAdaptativeGrid(TestCase):


    def plot_result(self, t, x, y):

        plt.figure()
        plt.plot(x[0], y[0,:],"k-")
        plt.plot(x[-1], y[-1,:],"bo-")
        plt.ylabel('y')
        plt.xlabel('distance')
        plt.title('Model1 Solution')
        plt.legend(["t={}".format(t[0]), "t={}".format(t[-1])])
        plt.show()

    def test_1(self):

        vL = 1.0
        vR = 0.5
        xi = 0.2
        N = 50

        lb = ModelParameter("lb", vL, "m/s")
        ub = ModelParameter("ub", vR, "m/s")

        t0 = np.linspace(0, 0.2, 2)
        tgrid = np.linspace(t0[0], t0[-1], 10)


        xa = np.linspace(0, 1., N)
        x0 = np.linspace(0, 1., N)

        uref = np.array([1])

        y0 = np.zeros(N)
        y0 = np.where(x0 < xi, vL, y0)
        y0 = np.where(x0 >= xi, vR, y0)

        NP = 50
        T = np.empty((len(tgrid)-1,))
        Y = np.empty((len(tgrid)-1, NP))
        X = np.empty((len(tgrid)-1, NP))

        x_adapt = x0
        y_adapt = y0

        for i in range(1,len(tgrid)):
            ti = tgrid[i-1]
            tf = tgrid[i]
            t = np.linspace(ti, tf, 4)
            print(ti, tf)

            x_adapt, y_adapt, y_linear = AdaptativeGrid.run(4, 12, NP, x_adapt, np.reshape(y_adapt, (N,1)), uref, eps = 1e-3, order = 5)
            model = Burgers(x_adapt, lb=lb, ub=ub, scheme=SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.CUBISTA2)
            t, y, yp = DasslSolver.run(model, t,  y_linear, display=False, rtol=1e-3, atol=1e-3)
            # time.sleep(1)
            print("..1.")

            y_adapt = y[-1,:]

            print("..2.")
            T[i-1] = t[-1]
            print("..3.")
            X[i-1,:] = x_adapt
            print("..4.")
            Y[i-1,:] = y_adapt
            print("..5.")

        time.sleep(2)
        self.plot_result(T, X, Y)

        self.assertTrue(True)
