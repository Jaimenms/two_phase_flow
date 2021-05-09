from unittest import TestCase
from models.two_phase_flow import TwoPhaseFlow
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.parameter import Parameters, ConstantParameter, TimeDependentParameter
from models.model.domain import Domain, Domains
from models.model.variable import Variable, Variables
from solvers.dassl_solver import DasslSolver


class TestTwoPhaseFlow(TestCase):

    def test_1(self, plot=True):

        N = 20
        t = np.linspace(0, 20, 10)
        x = np.linspace(0, 500.,N)

        alphaLi = 0.84
        alphaLf = 0.84
        Pi = 1200000
        Pf = 1000000
        qLi = 21
        qLf = 20
        qGi = 0.1
        qGf = 0.1

        D = ConstantParameter('D', 4, "in")
        epw = ConstantParameter('epw', 46, "um")
        g = ConstantParameter('g', 9.81, "m/s**2")
        rhoL = ConstantParameter('rhoL', 1000, "kg/m**3")
        muL = ConstantParameter('muL', 0.001, "Pa*s")
        drhoLdP = ConstantParameter('drhoLdP', 5e-7, "kg/m**3/Pa")
        rhoG = ConstantParameter('rhoG', 1, "kg/m**3")
        muG = ConstantParameter('muG', 0.00001, "Pa*s")
        drhoGdP = ConstantParameter('drhoGdP', 5e-4, "kg/m**3/Pa")
        z = ConstantParameter('z', np.zeros_like(x), "m")
        tspan = np.array([0, 1, 3, 10000])
        qL_lb = TimeDependentParameter('qL_lb', qLi * np.array([1, 1, 1.01, 1.01]), "kg/s", tspan=tspan)
        qG_lb = TimeDependentParameter('qG_lb', qGi * np.array([1, 1, 1.01, 1.01]), "kg/s", tspan=tspan)
        P_ub = TimeDependentParameter('P_ub', Pf * np.array([1, 1, 1, 1]), "Pa", tspan=tspan)

        parameters = Parameters((D, epw, g, rhoL, muL, drhoLdP, rhoG, muG, drhoGdP, z, qL_lb, qG_lb, P_ub))

        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        domains = Domains((x_domain,))

        qL0 = np.linspace(qLi, qLf, N)
        qL = Variable("qL", domains=(x_domain,), value=qL0, unit="kg/s")

        qG0 = np.linspace(qGi, qGf, N)
        qG = Variable("qG", domains=(x_domain,), value=qG0, unit="kg/s")

        alphaL0 = np.linspace(alphaLi, alphaLf, N)
        alphaL = Variable("alphaL", domains=(x_domain,), value=alphaL0, unit="")

        P0 = np.linspace(Pi, Pf, N)
        P = Variable("P", domains=(x_domain,), value=P0, unit="Pa")

        variables = Variables((qL, qG, alphaL, P))

        m = TwoPhaseFlow(
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
            self.assertAlmostEqual(y[-1,0], qLi)
