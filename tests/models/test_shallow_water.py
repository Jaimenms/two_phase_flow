from unittest import TestCase
from models.shallow_water import ShallowWater
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.parameter import Parameters, ConstantParameter, TimeDependentParameter
from models.model.domain import Domain, Domains
from models.model.variable import Variable, Variables
from solvers.dassl_solver import DasslSolver
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum


class TestShallowWater(TestCase):

    def test_1(self, plot=True):

        hL = 0.75
        hR = 0.25
        xi = 0.5
        N = 50

        t = np.linspace(0, 0.1, 4)

        parameters = Parameters(())

        x = np.linspace(0, 1., N)
        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")
        domains = Domains((x_domain,))

        v0 = np.zeros(N)
        v = Variable("v", domains=(x_domain,), value=v0, unit="kg/s")

        h0 = np.zeros(N)
        h0 = np.where(x < xi, hL, h0)
        h0 = np.where(x >= xi, hR, h0)
        h = Variable("h", domains=(x_domain,), value=h0, unit="Pa")

        variables = Variables((v, h))

        m = ShallowWater(
            domains=domains,
            variables=variables,
            parameters=parameters,
            scheme=SchemeM1FDMEnum.CENTRAL_N2, scheme_hrs=SchemeM1FDMEnum.CENTRAL_N2,
            flux_delimiter=FluxDelimiterEnum.SMART2
        )

        y0 = m.variables.read()

        t0, y0, yp0 = DasslSolver.run(m, None, y0, display=True, rtol=1e-3, atol=1e-6)

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], v0[0])