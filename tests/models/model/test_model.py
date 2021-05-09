from unittest import TestCase
import numpy as np
from models.model.variable import Variable
from models.model.variable import Domain, RegionEnum
from models.model.model import Model

class TestModel(TestCase):

    def test_1(self):
        x_val = np.linspace(0, 10, 11)
        x = Domain("x", value=x_val, unit="m", description="x coordinate")
        u = Variable("u-velocity", domains=(x,))
        value = np.linspace(0, 100, 11)
        value_ = Model.apply_regions(value, regions=(RegionEnum.OPEN_CLOSED,))
        self.assertEqual(len(value_), len(value)-1)
