from unittest import TestCase
import numpy as np
from methods.fdm.operations.second_gradient import SecondGradient
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum


class TestOperationSecondGradient(TestCase):

    def setUp(self):
        self.x4 = np.linspace(0,1,10)
        self.f4, self.d2f4 = self.x4 ** 4, 12 * self.x4 ** 2
        self.xs = np.linspace(0, 1., 100)
        self.fs = np.sin(2*3.1415*self.xs) + np.sin(3.1415*self.xs)/2
        self.dfs = 2*3.1415*np.cos(2*3.1415*self.xs) + 3.1415*np.cos(3.1415*self.xs)/2
        self.d2fs = -(2*3.1415)**2*np.sin(2*3.1415*self.xs) - 3.1415**2*np.sin(3.1415*self.xs)/2

    def subtests(self, df_correct, df_calc, places=7):
        for i, (df_correct_i, df_calc_i) in enumerate(zip(df_correct, df_calc)):
            with self.subTest(i=i):
                self.assertAlmostEqual(df_correct_i, df_calc_i, places=places)

    def test_central_n4(self):
        operator = SecondGradient(self.x4, scheme=SchemeM2FDMEnum.CENTRAL_N4)
        df_calc = operator(self.f4)
        self.subtests(self.d2f4, df_calc, places=7)

    def test_central_n4_2(self):
        operator = SecondGradient(self.xs, scheme=SchemeM2FDMEnum.CENTRAL_N4)
        df_calc = operator(self.fs)
        self.subtests(self.d2fs, df_calc, places=1)
