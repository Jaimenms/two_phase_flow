from unittest import TestCase
import numpy as np
from models.model.model_variable import ModelVariable
from models.model.model_variable import ModelDomain, ModelRegionEnum
from models.model.model import Model

class TestModel(TestCase):

    def test_1(self):
        x_val = np.linspace(0, 10, 11)
        x = ModelDomain("x", value=x_val, unit="m", description="x coordinate")
        u = ModelVariable("u-velocity", domains=(x,))
        value = np.linspace(0, 100, 11)
        value_ = Model.apply_regions(value,regions=(ModelRegionEnum.OPEN_CLOSED,))
        self.assertEqual(len(value_), len(value)-1)