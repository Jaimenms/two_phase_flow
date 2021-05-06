from unittest import TestCase
from solvers.dassl_solver import DasslSolver
from models.shallow_water import ShallowWater
import numpy as np
import matplotlib.pyplot as plt
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum


class TestShallowWater(TestCase):

    def case1(self):

        hL = 0.75
        hR = 0.25
        xi = 0.5
        N = 50

        t0 = np.linspace(0, 0.1, 4)
        x = np.linspace(0, 1., N)

        v0 = np.zeros(N)
        h0 = np.zeros(N)
        h0 = np.where(x < xi, hL, h0)
        h0 = np.where(x >= xi, hR, h0)

        model = ShallowWater(x, scheme=SchemeM1FDMEnum.CENTRAL_N2, scheme_hrs=SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.SMART2)

        y0 = np.concatenate([v0, h0], axis=None)

        return t0, y0, x, model


    def plot_result(self, t, x, y):

        fig, axs = plt.subplots(2)
        axs[0].plot(x, y[0,len(x):2*len(x)],"k-")
        axs[0].plot(x, y[-1,len(x):2*len(x)],"bo-")
        axs[0].set_ylabel('h')
        axs[0].set_xlabel('distance')
        axs[0].legend(["t={}".format(t[0]), "t={}".format(t[-1])])

        axs[1].plot(x, y[0,0*len(x):len(x)],"k-")
        axs[1].plot(x, y[-1,0*len(x):len(x)],"bo-")
        axs[1].set_ylabel('v')
        axs[1].set_xlabel('distance')
        axs[1].legend(["t={}".format(t[0]), "t={}".format(t[-1])])

        plt.show()


    def test_1(self, plot=True):

        t0, y0, x, model = self.case1()

        res, ires = model.residue(t0, y0, y0)

        t, y, yp = DasslSolver.run(model, t0, y0, display=True, rtol=1e-3, atol=1e-9)

        if plot:
            self.plot_result(t, x, y)
        else:
            self.assertAlmostEqual(t[-1], t0[-1])
            #self.assertAlmostEqual(y[-1][0], 0.135, places=3)
