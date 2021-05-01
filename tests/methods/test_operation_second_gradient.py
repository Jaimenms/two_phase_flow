from unittest import TestCase
import numpy as np
from methods.fdm.operations.second_gradient import SecondGradient
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum


class TestOperationSecondGradient(TestCase):

    def setUp(self):
        self.x = np.linspace(0,1,10)
        self.f4, self.d2f4 = self.x ** 4, 12 * self.x ** 2

    def subtests(self, df_correct, df_calc, places=7):
        for i, (df_correct_i, df_calc_i) in enumerate(zip(df_correct, df_calc)):
            with self.subTest(i=i):
                self.assertAlmostEqual(df_correct_i, df_calc_i, places=places)

    def test_central_n4(self):
        operator = SecondGradient(self.x, scheme=SchemeM2FDMEnum.CENTRAL_N4)
        df_calc = operator(self.f4)
        self.subtests(self.d2f4, df_calc, places=7)
