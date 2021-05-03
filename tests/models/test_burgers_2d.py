from unittest import TestCase
from models.burgers_2d import Burgers2D
import numpy as np


class TestBurgers2D(TestCase):

    def test_1(self):

        t = np.linspace(0, 5, 10)
        x1 = np.linspace(0, 2.,10)
        x2 = np.linspace(0, 2.,10)
        u = np.squeeze(np.zeros((10,10)))
        up = np.zeros_like(u)

        yp = np.concatenate((up, ), axis=None)
        y = np.concatenate((up, ), axis=None)

        m = Burgers2D(x1, x2)
        res, ires = m(t[0], y, yp)
        self.assertAlmostEqual(res[0], 0.0)
