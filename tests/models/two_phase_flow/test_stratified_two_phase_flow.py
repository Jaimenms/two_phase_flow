from unittest import TestCase
from models.two_phase_flow.two_phase_flow import TwoPhaseFlow
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.domain import Domain
from solvers.dassl_solver import DasslSolver


class TestStratifiedTwoPhaseFlow(TestCase):

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

        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        m = TwoPhaseFlow(
            x_domain=x_domain,
            flux_delimiter=FluxDelimiterEnum.SMART,
        )

        m.parameters['D'].set(4, "in")
        m.parameters['epw'].set(46, "um")
        m.parameters['g'].set(9.81, "m/s**2")
        m.parameters['rhoL'].set(1000, "kg/m**3")
        m.parameters['muL'].set(0.001, "Pa*s")
        m.parameters['drhoLdP'].set(5e-7, "kg/m**3/Pa")
        m.parameters['rhoG'].set(1, "kg/m**3")
        m.parameters['muG'].set(0.00001, "Pa*s")
        m.parameters['drhoGdP'].set(5e-4, "kg/m**3/Pa")
        m.parameters['tetha'].set(np.zeros_like(x), "")
        m.parameters['qL_lb'].set(value_span=qLi * np.array([1, 1, 1.01, 1.01]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 1, 3, 10000]))
        m.parameters['qG_lb'].set(value_span=qGi * np.array([1, 1, 1.01, 1.01]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 1, 3, 10000]))
        m.parameters['P_ub'].set(value_span=Pf * np.array([1, 1, 1, 1]), value_unit="Pa", time_unit="s", time_span=np.array([0, 1, 3, 10000]))

        qL0 = np.linspace(qLi, qLf, N)
        m.variables["qL"].set_ic(value=qL0, unit="kg/s")
        qG0 = np.linspace(qGi, qGf, N)
        m.variables["qG"].set_ic(value=qG0, unit="kg/s")
        alphaL0 = np.linspace(alphaLi, alphaLf, N)
        m.variables["alphaL"].set_ic(value=alphaL0, unit="")
        P0 = np.linspace(Pi, Pf, N)
        m.variables["P"].set_ic(value=P0, unit="Pa")

        y0 = m.variables.get_ic_array()

        t0, y0, yp0 = DasslSolver.run(m, None, y0, display=True, rtol=1e-3, atol=1e-6)

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], qLi)
