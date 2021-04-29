from unittest import TestCase
from models.model.model_domain import ModelDomain
import numpy as np


class TestModelDomain(TestCase):

    def test_case1(self):
        value = np.linspace(0, 10, 11)
        x = ModelDomain("x", value=np.linspace(0,10,11), unit="cm", description="x coordinate")

        self.assertEqual(x.base_unit, "meter")
        self.assertEqual(value[-1], x[-1]*100)
        self.assertEqual(len(x), len(value))
