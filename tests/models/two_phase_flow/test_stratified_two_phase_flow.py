from unittest import TestCase
from models.two_phase_flow.two_phase_flow import TwoPhaseFlow
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.domain import Domain
from solvers.dassl_solver import DasslSolver


class TestStratifiedTwoPhaseFlow(TestCase):

    def test_1(self, plot=True):

        N = 25
        t = np.linspace(0, 100, 11)
        x = np.linspace(0, 100.,N)

        Pr=1e5

        alphaLi = 0.46
        Pi = 520000/Pr
        Pf = 500000/Pr
        qLi = 11
        qGi = 0.5


        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        m = TwoPhaseFlow(
            x_domain=x_domain,
            flux_delimiter=FluxDelimiterEnum.MINMOD2,
        )

        m.parameters['Pr'].set(Pr, "Pa")
        m.parameters['D'].set(8, "in")
        m.parameters['epw'].set(46, "um")
        m.parameters['g'].set(9.81, "m/s**2")
        m.parameters['rhoL'].set(1000, "kg/m**3")
        m.parameters['muL'].set(0.001, "Pa*s")
        m.parameters['drhoLdP'].set(5e-7, "kg/m**3/Pa")
        m.parameters['rhoG'].set(1, "kg/m**3")
        m.parameters['muG'].set(0.00001, "Pa*s")
        m.parameters['drhoGdP'].set(5e-4, "kg/m**3/Pa")
        m.parameters['tetha'].set(np.zeros_like(x), "")
        m.parameters['qL_lb'].set(value_span=qLi * np.array([1, 1, 1.1, 1.1]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 20, 60, 10000]))
        m.parameters['qG_lb'].set(value_span=qGi * np.array([1, 1, 1., 1.]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 2, 3, 10000]))
        m.parameters['P_ub'].set(value_span=Pf * Pr * np.array([1, 1, 1, 1]), value_unit="Pa", time_unit="s", time_span=np.array([0, 1, 3, 10000]))

        qL0 = np.linspace(qLi, qLi, N)
        m.variables["qL"].set_ic(value=qL0, unit="kg/s")
        qG0 = np.linspace(qGi, qGi, N)
        m.variables["qG"].set_ic(value=qG0, unit="kg/s")
        alphaL0 = np.linspace(alphaLi, alphaLi, N)
        m.variables["alphaL"].set_ic(value=alphaL0, unit="")
        P0 = np.linspace(Pi, Pf, N)
        m.variables["P"].set_ic(value=P0, unit="Pa/Pa")

        y0 = m.variables.get_ic_array()

        out = m.jacobian(0, y0, np.zeros_like(y0), 1.0, None)

        t0, y0, yp0 = DasslSolver.run(m, None, y0, display=True, rtol=1e-3, atol=1e-6)

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6, user_jacobian=True)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], qLi)


    def test_2(self, plot=True):

        N = 20
        t = np.linspace(0, 300, 11)
        x = np.linspace(0, 100.,N)

        Pr=1e5

        alphaLi = 0.46
        Pi = 520000/Pr
        Pf = 500000/Pr
        qLi = 11
        qGi = 0.5

        x_domain = Domain("x", value=x, unit="m", description="x1 coordinate")

        m = TwoPhaseFlow(
            x_domain=x_domain,
            flux_delimiter=FluxDelimiterEnum.SMART,
        )


        tetha = np.ones_like(x)
        tetha[x<50] = 2*np.pi*(-0.1/360)
        tetha[x>=50] = 2*np.pi*(1/360)
        m.parameters['tetha'].set(tetha, "")

        m.parameters['Pr'].set(Pr, "Pa")
        m.parameters['D'].set(8, "in")
        m.parameters['epw'].set(46, "um")
        m.parameters['g'].set(9.81, "m/s**2")
        m.parameters['rhoL'].set(1000, "kg/m**3")
        m.parameters['muL'].set(0.001, "Pa*s")
        m.parameters['drhoLdP'].set(5e-7, "kg/m**3/Pa")
        m.parameters['rhoG'].set(1, "kg/m**3")
        m.parameters['muG'].set(0.00001, "Pa*s")
        m.parameters['drhoGdP'].set(5e-4, "kg/m**3/Pa")
        m.parameters['qL_lb'].set(value_span=qLi * np.array([1, 1, 1.1, 1.1]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 20, 60, 10000]))
        m.parameters['qG_lb'].set(value_span=qGi * np.array([1, 1, 1., 1.]), value_unit="kg/s", time_unit="s", time_span=np.array([0, 2, 3, 10000]))
        m.parameters['P_ub'].set(value_span=Pf * Pr * np.array([1, 1, 1, 1]), value_unit="Pa", time_unit="s", time_span=np.array([0, 1, 3, 10000]))

        qL0 = np.linspace(qLi, qLi, N)
        m.variables["qL"].set_ic(value=qL0, unit="kg/s")
        qG0 = np.linspace(qGi, qGi, N)
        m.variables["qG"].set_ic(value=qG0, unit="kg/s")
        alphaL0 = np.linspace(alphaLi, alphaLi, N)
        m.variables["alphaL"].set_ic(value=alphaL0, unit="")
        P0 = np.linspace(Pi, Pf, N)
        m.variables["P"].set_ic(value=P0, unit="Pa/Pa")

        y0 = m.variables.get_ic_array()

        out = m.jacobian(0, y0, np.zeros_like(y0), 1.0, None)

        t0, y0, yp0 = DasslSolver.run(m, None, y0, display=True, rtol=1e-3, atol=1e-4)

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-4, user_jacobian=True)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], qLi)
