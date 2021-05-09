from unittest import TestCase
from models.burgers import Burgers
import numpy as np
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from models.model.parameter import Parameter


class TestBurgers(TestCase):

    def test_1(self):

        t = np.linspace(0, 5, 10)
        x = np.linspace(0, 2.,100)
        y = np.zeros(100)
        yp = np.zeros(100)

        lb = Parameter("lb", 0.0, "m/s")
        ub = Parameter("ub", 0.0, "m/s")

        m = Burgers(x, lb=lb, ub=ub, flux_delimiter=FluxDelimiterEnum.SMART)
        res, ires = m(t[0], y, yp)
        self.assertAlmostEqual(res[0], 0.0)
