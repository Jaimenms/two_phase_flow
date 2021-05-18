from unittest import TestCase
from models.shallow_water import ShallowWater
import numpy as np
from models.model.domain import Domain
from solvers.dassl_solver import DasslSolver
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum


class TestShallowWater(TestCase):

    def test_jac(self):

        hL = 0.75
        hR = 0.25
        xi = 0.5
        N = 8
        t = np.linspace(0, 0.1, 4)

        x = np.linspace(0, 1., N)
        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        m = ShallowWater(
            x_domain=x_domain,
            scheme=SchemeM1FDMEnum.CENTRAL_N2, scheme_hrs=SchemeM1FDMEnum.CENTRAL_N2,
            flux_delimiter=FluxDelimiterEnum.SMART2
        )

        v0 = np.zeros(N)
        m.variables["v"].set_ic(v0,"m/s")

        h0 = np.zeros(N)
        h0 = np.where(x < xi, hL, h0)
        h0 = np.where(x >= xi, hR, h0)
        m.variables["h"].set_ic(h0,"m")

        y0 = m.variables.get_ic_array()

        ou1, ou2 = m.jacobian(0, y0, np.zeros_like(y0), 1.0,None)

        self.assertTrue(True)


    def test_1(self, plot=True):

        # TODO  - This is not working

        hL = 0.75
        hR = 0.25
        xi = 0.5
        N = 200
        t = np.linspace(0, 0.1, 4)

        x = np.linspace(0, 1., N)
        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        m = ShallowWater(
            x_domain=x_domain,
            scheme=SchemeM1FDMEnum.CENTRAL_N4, scheme_hrs=SchemeM1FDMEnum.CENTRAL_N4,
            flux_delimiter=FluxDelimiterEnum.MINMOD2
        )

        v0 = np.zeros(N)
        m.variables["v"].set_ic(v0,"m/s")

        h0 = np.zeros(N)
        h0 = np.where(x < xi, hL, h0)
        h0 = np.where(x >= xi, hR, h0)
        m.variables["h"].set_ic(h0,"m")

        y0 = m.variables.get_ic_array()

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-3, user_jacobian=True)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], v0[0])