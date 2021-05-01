from unittest import TestCase
import numpy as np
from methods.fdm.fdm_error import FDMError
from methods.fdm.operations.first_gradient import FirstGradient
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum


class TestOperationFirstGradient(TestCase):

    def setUp(self):
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
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N2)
        df_calc = operator(self.f2)
        self.subtests(self.df2, df_calc, places=7)

    def test_central_n4(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N4)
        df_calc = operator(self.f4)
        self.subtests(self.df4, df_calc, places=7)

    def test_central_n6(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N6)
        df_calc = operator(self.f6)
        self.subtests(self.df6, df_calc, places=7)

    def test_central_n8(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.CENTRAL_N8)
        df_calc = operator(self.f8)
        self.subtests(self.df8, df_calc, places=7)

    def test_central_n8_error(self):
        x = np.logspace(0, 1, 5)
        with self.assertRaises(FDMError) as context:
            operator = FirstGradient(x, scheme=SchemeM1FDMEnum.CENTRAL_N8)

    def test_upwind_n2(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N2)
        df_calc = operator(self.f2)
        self.subtests(self.df2, df_calc, places=7)

    def test_upwind_n4(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N4)
        df_calc = operator(self.f4)
        self.subtests(self.df4, df_calc, places=7)

    def test_upwind_n6(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N6)
        df_calc = operator(self.f6)
        self.subtests(self.df6, df_calc, places=7)

    def test_upwind_n8(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.UPWIND_N8)
        df_calc = operator(self.f8)
        self.subtests(self.df8, df_calc, places=7)

    def test_downwind_n2(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N2)
        df_calc = operator(self.f2)
        self.subtests(self.df2, df_calc, places=7)

    def test_downwind_n4(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N4)
        df_calc = operator(self.f4)
        self.subtests(self.df4, df_calc, places=7)

    def test_downwind_n6(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N6)
        df_calc = operator(self.f6)
        self.subtests(self.df6, df_calc, places=7)

    def test_downwind_n8(self):
        operator = FirstGradient(self.x, scheme=SchemeM1FDMEnum.DOWNWIND_N8)
        df_calc = operator(self.f8)
        self.subtests(self.df8, df_calc, places=7)

    def test_central_n2_2d(self):

        x = np.linspace(0, 5, 6)
        y = np.linspace(0, 7, 8)
        X, Y = np.meshgrid(x,y, indexing='ij')

        f1 = X**2
        f2 = Y**3

        df1dx = 2*X
        df1dy = np.zeros_like(X)
        df2dx = np.zeros_like(Y)
        df2dy = 3*Y**2

        grad_x = FirstGradient(x, axis=0, scheme=SchemeM1FDMEnum.CENTRAL_N4)
        grad_y = FirstGradient(y, axis=1, scheme=SchemeM1FDMEnum.CENTRAL_N4)

        df1dx_calc = grad_x(f1)
        df1dy_calc = grad_y(f1)
        df2dx_calc = grad_x(f2)
        df2dy_calc = grad_y(f2)

        np.testing.assert_array_almost_equal(df1dx, df1dx_calc, decimal=6)
        np.testing.assert_array_almost_equal(df1dy, df1dy_calc, decimal=6)
        np.testing.assert_array_almost_equal(df2dx, df2dx_calc, decimal=6)
        np.testing.assert_array_almost_equal(df2dy, df2dy_calc, decimal=6)


    def test_central_n2_3d(self):

        x = np.linspace(0, 5, 6)
        y = np.linspace(10, 17, 8)
        z = np.linspace(20, 30, 11)
        X, Y, Z = np.meshgrid(x,y,z, indexing='ij')

        f1 = X**2
        f2 = Y**3
        f3 = Z**4

        df1dx = 2*X
        df1dy = np.zeros_like(X)
        df1dz = np.zeros_like(X)
        df2dx = np.zeros_like(Y)
        df2dy = 3*Y**2
        df2dz = np.zeros_like(Y)
        df3dx = np.zeros_like(Y)
        df3dy = np.zeros_like(Y)
        df3dz = 4*Z**3

        grad_x = FirstGradient(x, axis=0, scheme=SchemeM1FDMEnum.CENTRAL_N4)
        grad_y = FirstGradient(y, axis=1, scheme=SchemeM1FDMEnum.CENTRAL_N4)
        grad_z = FirstGradient(z, axis=2, scheme=SchemeM1FDMEnum.CENTRAL_N4)

        df1dx_calc = grad_x(f1)
        df1dy_calc = grad_y(f1)
        df1dz_calc = grad_z(f1)
        df2dx_calc = grad_x(f2)
        df2dy_calc = grad_y(f2)
        df2dz_calc = grad_z(f2)
        df3dx_calc = grad_x(f3)
        df3dy_calc = grad_y(f3)
        df3dz_calc = grad_z(f3)

        np.testing.assert_array_almost_equal(df1dx, df1dx_calc, decimal=6)
        np.testing.assert_array_almost_equal(df1dy, df1dy_calc, decimal=6)
        np.testing.assert_array_almost_equal(df1dz, df1dz_calc, decimal=6)
        np.testing.assert_array_almost_equal(df2dx, df2dx_calc, decimal=6)
        np.testing.assert_array_almost_equal(df2dy, df2dy_calc, decimal=6)
        np.testing.assert_array_almost_equal(df2dz, df2dz_calc, decimal=6)
        np.testing.assert_array_almost_equal(df3dx, df3dx_calc, decimal=6)
        np.testing.assert_array_almost_equal(df3dy, df3dy_calc, decimal=6)
        np.testing.assert_array_almost_equal(df3dz, df3dz_calc, decimal=6)
