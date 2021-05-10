from unittest import TestCase
from models.single_phase_flow import SinglePhaseFlow
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.domain import Domain
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

        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        m = SinglePhaseFlow(
            x_domain=x_domain,
            flux_delimiter=FluxDelimiterEnum.SMART,
        )

        m.parameters["D"].set(4, "in")
        m.parameters['g'].set(9.81, "m/s**2")
        m.parameters['epw'].set(46, "um")
        m.parameters['mu'].set(0.001, "Pa*s")
        m.parameters['rho'].set(1000, "kg/m**3")
        m.parameters['drhodP'].set(5e-7, "kg/m**3/Pa")
        m.parameters['z'].set(np.zeros_like(x), "m")
        m.parameters['q_lb'].set(value_span=qi * np.array([1, 1, 1.01, 1.01]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 1, 3, 10000]))
        m.parameters['P_ub'].set(value_span=Pf * np.array([1, 1, 1, 1]), value_unit="Pa", time_unit="s", time_span=np.array([0, 1, 3, 10000]))

        q0 = np.linspace(qi, qf, N)
        m.variables["q"].set_ic(value=q0, unit="kg/s")

        P0 = np.linspace(Pi, Pf, N)
        m.variables["P"].set_ic(value=P0, unit="Pa")

        y0 = m.variables.get_ic_array()

        t0, y0, yp0 = DasslSolver.run(m, None, y0, display=True, rtol=1e-3, atol=1e-6)

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], qi)