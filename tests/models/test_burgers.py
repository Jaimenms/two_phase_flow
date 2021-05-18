from unittest import TestCase
from models.burgers import Burgers
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.domain import Domain
from solvers.dassl_solver import DasslSolver
from scipy.interpolate import interp1d
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum


class TestBurgers(TestCase):

    def test_1(self, plot=True):

        # TODO - Does not converge for N > 200

        N = 100
        t = np.linspace(0, 0.8, 110)
        tresp = np.linspace(0, 0.8, 4)
        x = np.linspace(0, 1., N)
        u0 = np.sin(2*3.1415*x) + np.sin(3.1415*x)/2
        deltax = x[1] - x[0]
        vmax = np.max(u0)
        Cmax = 1
        deltatmax = Cmax*deltax/vmax

        m = Burgers(
            x_domain = Domain("x", x, "m"),
            flux_delimiter=FluxDelimiterEnum.SMART2,
            scheme=SchemeM1FDMEnum.CENTRAL_N2, scheme_hrs=SchemeM1FDMEnum.CENTRAL_N6,
        )

        m.parameters["visc"].set(0.001, "m**2/s")
        m.parameters["lb"].set(None, "m/s")
        m.parameters["ub"].set(None, "m/s")
        #m.jacobian = None

        m.variables["u"].set_ic(u0,"m/s")
        y0 = m.variables.get_ic_array()

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6, user_jacobian=True)

        fresp = interp1d(t, y, axis=0)
        yresp = fresp(tresp)

        if plot:
            m.plot_result(tresp, yresp)
        else:
            self.assertAlmostEqual(y[-1,0], u0[0])


    def test_2(self, plot=False):

        N = 1000
        #t = np.linspace(0, 0.6, 600)
        t = np.linspace(0, 0.6, 4)

        vL = 1.0
        vR = 0.5
        xi = 0.2
        x = np.linspace(0, 1., N)
        deltax = x[1] - x[0]
        vmax = vL
        Cmax = 1
        deltatmax = Cmax*deltax/vmax

        m = Burgers(
            x_domain = Domain("x", x, "m"),
            flux_delimiter=FluxDelimiterEnum.SMART2
        )

        m.parameters["visc"].set(0.001, "m**2/s")
        m.parameters["lb"].set(vL, "m/s")
        m.parameters["ub"].set(vR, "m/s")

        u0 = np.zeros(N)
        u0 = np.where(x < xi, vL, u0)
        u0 = np.where(x >= xi, vR, u0)
        m.variables["u"].set_ic(u0,"m/s")

        y0 = m.variables.get_ic_array()
        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6, user_jacobian=True)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], u0[0])
