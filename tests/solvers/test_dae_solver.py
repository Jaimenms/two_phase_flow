from unittest import TestCase
from solvers.dae_solver import DAESolver
from models.burgers import Burgers
import numpy as np
import matplotlib.pyplot as plt
from methods.fdm.fdm_mixin import FluxDelimiterEnum, FDMEnum

class TestDAESolver(TestCase):

    def test_burgers(self):

        vL = 1.0
        vR = 0.5
        xi = 0.2
        N = 200

        t0 = np.linspace(0, 0.5, 2)
        x = np.linspace(0, 1., N)

        KIND_IC = "sine"
        if KIND_IC == "step":
            y0 = np.zeros(N)
            y0 = np.where(x < xi, vL, y0)
            y0 = np.where(x >= xi, vR, y0)
        else:
            y0 = np.sin(2*3.1415*x) + np.sin(3.1415*x)/2

        m = Burgers(x, order=FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.MINMOD)

        it1=0
        it2=-1
        t, y, yp = DAESolver().run(m.residue, t0, y0, None, None, 1e-6, 1e-8)
        self.assertAlmostEqual(t[-1], t0[-1])
        #self.assertAlmostEqual(y[-1][0], 0.135, places=3)
        print(y)
        plt.figure(2)
        plt.plot(x, y[it1,:],"k-")
        plt.plot(x, y[it2,:],"b-")
        plt.ylabel('y')
        plt.xlabel('distance')
        plt.title('Model1 Solution')
        plt.legend(["t={}".format(t[it1]), "t={}".format(t[it2])])
        plt.show()
