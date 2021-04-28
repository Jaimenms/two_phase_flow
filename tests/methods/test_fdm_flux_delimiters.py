from unittest import TestCase
import numpy as np
from methods.fdm.fdm_mixin import FDMMixin, FDMEnum, FluxDelimiterEnum
from methods.fdm.fdm_error import FDMError


class TestFdmFluxDelimiters(TestCase):

    def setUp(self):
        self.m = FDMMixin()
        self.x = np.linspace(0,1,10)
        self.f0, self.df0 = self.x ** 0, 0 * self.x ** 0
        self.f1, self.df1 = self.x ** 1, 1 * self.x ** 0
        self.f2, self.df2 = self.x ** 2, 2 * self.x ** 1
        self.f3, self.df3 = self.x ** 3, 3 * self.x ** 2
        self.f4, self.df4 = self.x ** 4, 4 * self.x ** 3
        self.f6, self.df6 = self.x ** 6, 6 * self.x ** 5
        self.f8, self.df8 = self.x ** 8, 8 * self.x ** 7
        self.fe, self.dfe = np.exp(self.x/100), 1/100 * np.exp(self.x/100)
        self.f2a, self.df2a = (self.x-1) ** 2, 2 * (self.x-1) ** 1

        self.fh = (1-np.heaviside(np.linspace(-10,10,10), 4))*0.5 + 0.5
        self.dfh = np.zeros_like(self.fh)

    def subtests(self, df_correct, df_calc, places=7):
        for i, (df_correct_i, df_calc_i) in enumerate(zip(df_correct, df_calc)):
            with self.subTest(i=i):
                self.assertAlmostEqual(df_correct_i, df_calc_i, places=places)

    def test_hrs_1(self):
        operator = self.m.Gradient(self.x, order=FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.MINMOD)
        df_calc = operator(self.f1, np.ones_like(self.f2))
        self.subtests(self.df1, df_calc, places=7)

    def test_hrs_2(self):
        operator = self.m.Gradient(self.x,order=FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.SMART)
        df_calc = operator(self.fh, np.ones_like(self.fh))
        self.subtests(self.dfh[0:5], df_calc[0:5], places=7)
        self.subtests(self.dfh[6:], df_calc[6:], places=7)
