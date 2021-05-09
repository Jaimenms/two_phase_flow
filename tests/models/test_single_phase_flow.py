from unittest import TestCase
from models.single_phase_flow import SinglePhaseFlow
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.parameter import Parameters, ConstantParameter, TimeDependentParameter
from models.model.domain import Domain, Domains
from models.model.variable import Variable, Variables
from solvers.dassl_solver import DasslSolver


class TestSinglePhaseFlow(TestCase):

    def test_1(self, plot=True):

        N = 200
        t = np.linspace(0, 20, 10)
        x = np.linspace(0, 500.,N)

        Pi = 1200000
        Pf = 1000000
        qi = 21
        qf = 20

        D = ConstantParameter('D', 4, "in")
        g = ConstantParameter('g', 9.81, "m/s**2")
        epw = ConstantParameter('epw', 46, "um")
        mu = ConstantParameter('mu', 0.001, "Pa*s")
        rho = ConstantParameter('rho', 1000, "kg/m**3")
        drhodP = ConstantParameter('drhodP', 5e-7, "kg/m**3/Pa")
        z = ConstantParameter('z', np.zeros_like(x), "m")

        tspan = np.array([0, 1, 3, 10000])
        q_lb = TimeDependentParameter('q_lb', qi * np.array([1, 1, 1.01, 1.01]), "kg/s", tspan=tspan)
        P_ub = TimeDependentParameter('P_ub', Pf * np.array([1, 1, 1, 1]), "Pa", tspan=tspan)

        parameters = Parameters((D, epw, g, rho, mu, drhodP, rho, mu, z, q_lb, P_ub))
        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")
        domains = Domains((x_domain,))

        q0 = np.linspace(qi, qf, N)
        q = Variable("q", domains=(x_domain,), value=q0, unit="kg/s")

        P0 = np.linspace(Pi, Pf, N)
        P = Variable("P", domains=(x_domain,), value=P0, unit="Pa")

        variables = Variables((q, P))

        m = SinglePhaseFlow(
            domains=domains,
            variables=variables,
            parameters=parameters,
            flux_delimiter=FluxDelimiterEnum.SMART,
        )

        y0 = m.variables.read()

        t0, y0, yp0 = DasslSolver.run(m, None, y0, display=True, rtol=1e-3, atol=1e-6)

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], qi)