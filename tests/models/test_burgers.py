from unittest import TestCase
from models.burgers import Burgers
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.domain import Domain
from solvers.dassl_solver import DasslSolver


class TestBurgers(TestCase):

    def test_1(self, plot=True):

        N = 100
        t = np.linspace(0, 0.5, 4)
        x = np.linspace(0, 1., N)
        u0 = np.sin(2*3.1415*x) + np.sin(3.1415*x)/2

        m = Burgers(
            x_domain = Domain("x", x, "m"),
            flux_delimiter=FluxDelimiterEnum.SMART
        )

        m.parameters["visc"].set(0.001, "m**2/s")
        m.variables["u"].set_ic(u0,"m/s")
        y0 = m.variables.get_ic_array()

        m.residue(t,y0,np.zeros_like(y0))

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], u0[0])


    def test_2(self, plot=True):

        N = 100
        t = np.linspace(0, 0.6, 4)

        vL = 1.0
        vR = 0.5
        xi = 0.2
        x = np.linspace(0, 1., N)

        m = Burgers(
            x_domain = Domain("x", x, "m"),
            flux_delimiter=FluxDelimiterEnum.SMART
        )

        m.parameters["visc"].set(0.001, "m**2/s")
        m.parameters["lb"].set(vL, "m/s")
        m.parameters["ub"].set(vR, "m/s")

        u0 = np.zeros(N)
        u0 = np.where(x < xi, vL, u0)
        u0 = np.where(x >= xi, vR, u0)
        m.variables["u"].set_ic(u0,"m/s")

        y0 = m.variables.get_ic_array()

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], u0[0])
