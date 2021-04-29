from unittest import TestCase
import numpy as np
from methods.fdm.fdm_mixin import FDMMixin, SchemeM1FDMEnum, FluxDelimiterEnum
from methods.fdm.fdm_error import FDMError


class TestFdmScuemes(TestCase):

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

    def test_central_n2(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N2)
        df_calc = operator(self.f2)
        self.subtests(self.df2, df_calc, places=7)

    def test_central_n4(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N4)
        df_calc = operator(self.f4)
        self.subtests(self.df4, df_calc, places=7)

    def test_central_n6(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N6)
        df_calc = operator(self.f6)
        self.subtests(self.df6, df_calc, places=7)

    def test_central_n8(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N8)
        df_calc = operator(self.f8)
        self.subtests(self.df8, df_calc, places=7)

    def test_central_n8_error(self):
        x = np.logspace(0, 1, 5)
        with self.assertRaises(FDMError) as context:
            operator = self.m.Gradient(x, scheme=SchemeM1FDMEnum.CENTRAL_N8)

    def test_upwind_n2(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N2)
        df_calc = operator(self.f2)
        self.subtests(self.df2, df_calc, places=7)

    def test_upwind_n4(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N4)
        df_calc = operator(self.f4)
        self.subtests(self.df4, df_calc, places=7)

    def test_upwind_n6(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N6)
        df_calc = operator(self.f6)
        self.subtests(self.df6, df_calc, places=7)

    def test_upwind_n8(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N8)
        df_calc = operator(self.f8)
        self.subtests(self.df8, df_calc, places=7)

    def test_downwind_n2(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N2)
        df_calc = operator(self.f2)
        self.subtests(self.df2, df_calc, places=7)

    def test_downwind_n4(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N4)
        df_calc = operator(self.f4)
        self.subtests(self.df4, df_calc, places=7)

    def test_downwind_n6(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N6)
        df_calc = operator(self.f6)
        self.subtests(self.df6, df_calc, places=7)

    def test_downwind_n8(self):
        operator = self.m.Gradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N8)
        df_calc = operator(self.f8)
        self.subtests(self.df8, df_calc, places=7)