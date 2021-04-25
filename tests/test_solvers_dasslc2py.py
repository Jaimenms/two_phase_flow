from unittest import TestCase
import dasslc
import numpy as np

class TestDasslc2py(TestCase):

    def test_model0(self):

        def model0(t, y, yp):
            res = np.empty(1)
            res[0] = yp[0] + 2 * y[0]
            ires = 0
            return res, ires

        t0 = np.array([1, 2])
        y0 = np.array([1])
        t, y, yp = dasslc.solve(model0, t0, y0)
        self.assertAlmostEqual(t[-1], 2, places=1)
        self.assertAlmostEqual(y[-1][0], 0.135, places=3)
