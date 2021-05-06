from unittest import TestCase
from solvers.dassl_solver import DasslSolver
from models.burgers_2d import Burgers2D
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum

class TestDasslSolverBurgers2D(TestCase):

    def case(self):

        N = 25

        #t0 = np.linspace(0, 0.158, 2)
        t0 = np.linspace(0, 0.5, 10)
        x1 = np.linspace(-1, 1., N)
        x2 = np.linspace(-1, 1., N)

        X1, X2 = np.meshgrid(x1, x2, indexing="ij")

        U = np.zeros_like(X1)
        U[(X1 + 0.2)**2 + (X2 + 0.2)**2 <=0.4] = 1.0

        y0 = np.concatenate((U, ), axis=None)
        model = Burgers2D(x1, x2, scheme=SchemeM1FDMEnum.CENTRAL_N4, flux_delimiter=FluxDelimiterEnum.SMART2)

        return t0, y0, x1, x2, X1, X2, model

    def plot_result(self, t, x1, x2, y):

        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        Y = np.reshape(y, (len(t), len(x1), len(x2)))

        ax.plot_surface(X1, X2, Y[-1])
        plt.ylabel('x2')
        plt.xlabel('x1')
        plt.title('Model1 Solution')
        plt.show()

    def test_burgers(self, plot=True):

        t0, y0, x1, x2, X1, X2, model = self.case()

        t, y, yp = DasslSolver.run(model, t0, y0, display=True, rtol=1e-3, atol=1e-3)

        if plot:
            self.plot_result(t, x1, x2, y)
        else:
            self.assertAlmostEqual(t[-1], t0[-1])
            #self.assertAlmostEqual(y[-1][0], 0.135, places=3)
