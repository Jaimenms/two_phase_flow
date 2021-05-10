from unittest import TestCase
from models.burgers_2d import Burgers2D
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.domain import Domain
from solvers.dassl_solver import DasslSolver


class TestBurgers2D(TestCase):

    def test_1(self, plot = True):

        N = 10
        t = np.linspace(0, 0.5, 4)
        x1 = np.linspace(-1, 1., N)
        x2 = np.linspace(-1, 1., N)

        x1_domain = Domain("x1", value=x1, unit="m", description="x1 coordinate")
        x2_domain = Domain("x2", value=x2, unit="m", description="x2 coordinate")

        m = Burgers2D(
            x1_domain=x1_domain,
            x2_domain=x2_domain,
            flux_delimiter=FluxDelimiterEnum.SMART
        )

        m.parameters["visc"].set(0.001, "m**2/s")

        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        u0 = np.zeros_like(X1)
        u0[(X1 + 0.2)**2 + (X2 + 0.2)**2 <=0.4] = 1.0
        m.variables["u"].set_ic(u0,"m/s")

        y0 = m.variables.get_ic_array()

        t, y, yp = DasslSolver.run(m, t, y0, display=True, rtol=1e-3, atol=1e-6)

        if plot:
            m.plot_result(t, y)
        else:
            self.assertAlmostEqual(y[-1,0], u0[0,0])

