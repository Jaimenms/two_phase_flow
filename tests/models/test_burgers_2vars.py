from unittest import TestCase
from models.burgers_2vars import Burgers
import numpy as np


class TestBurgers(TestCase):

    def test_1(self):

        t = np.linspace(0, 5, 10)
        x1 = np.linspace(0, 2.,10)
        x2 = np.linspace(0, 2.,10)
        y = np.squeeze(np.zeros((10,10)))
        yp = np.squeeze(np.zeros((10,10)))
        m = Burgers(x1, x2)
        res, ires = m(t[0], y, yp)
        self.assertAlmostEqual(res[0], 0.0)
