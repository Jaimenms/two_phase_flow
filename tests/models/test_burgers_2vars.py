from unittest import TestCase
from models.burgers_2vars import Burgers
import numpy as np


class TestBurgers(TestCase):

    def test_1(self):

        t = np.linspace(0, 5, 10)
        x1 = np.linspace(0, 2.,10)
        x2 = np.linspace(0, 2.,10)
        u = np.squeeze(np.zeros((10,10)))
        v = np.squeeze(np.ones((10,10)))
        up = np.zeros_like(u)
        vp = np.zeros_like(v)

        yp = np.concatenate((up, vp), axis=None)
        y = np.concatenate((up, vp), axis=None)

        m = Burgers(x1, x2)
        res, ires = m(t[0], y, yp)
        self.assertAlmostEqual(res[0], 0.0)
