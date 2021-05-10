from unittest import TestCase
from models.burgers import Burgers
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.parameter import Parameters, ConstantParameter, TimeDependentParameter
from models.model.domain import Domain, Domains
from models.model.variable import Variable, Variables
from solvers.dassl_solver import DasslSolver


class TestBurgers(TestCase):

    def test_1(self, plot=True):

        N = 100
        t = np.linspace(0, 0.5, 4)
        x = np.linspace(0, 1., N)

        lb = ConstantParameter('lb', None, "m/s")
        ub = ConstantParameter('ub', None, "m/s")

        parameters = Parameters((lb, ub))

        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")
        domains = Domains((x_domain,))

        u0 = np.sin(2*3.1415*x) + np.sin(3.1415*x)/2
        u = Variable("u", domains=(x_domain,), value=u0, unit="kg/s")

        variables = Variables((u,))

        m = Burgers(
            domains=domains,
            variables=variables,
            parameters=parameters,
            flux_delimiter=FluxDelimiterEnum.SMART
        )

        y0 = m.variables.read()

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

        lb = ConstantParameter('lb', vL, "m/s")
        ub = ConstantParameter('ub', vR, "m/s")

        parameters = Parameters((lb, ub))

        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")
        domains = Domains((x_domain,))

        u0 = np.zeros(N)
        u0 = np.where(x < xi, vL, u0)
        u0 = np.where(x >= xi, vR, u0)
        u = Variable("u", domains=(x_domain,), value=u0, unit="kg/s")

        variables = Variables((u,))

        m = Burgers(
            domains=domains,
            variables=variables,
            parameters=parameters,
            flux_delimiter=FluxDelimiterEnum.SMART
        )

        y0 = m.variables.read()

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], u0[0])