from unittest import TestCase
from models.burgers import Burgers
import numpy as np


class TestBurgers(TestCase):

    def test_1(self):

        t = np.linspace(0, 5, 10)
        x = np.linspace(0, 2.,100)
        y = np.zeros(100)
        yp = np.zeros(100)
        m = Burgers(x)
        res, ires = m(t[0], y, yp)
        self.assertAlmostEqual(res[0], 0.0)
