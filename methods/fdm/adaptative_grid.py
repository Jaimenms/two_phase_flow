import numpy as np
from scipy import interpolate


class AdaptativeGrid:

    @staticmethod
    def run(M: int, N: int, NPmax: int, x0: np.ndarray, u0: np.ndarray, uref: np.ndarray, eps: float = 1e-1, order: int = 5):
        """

        :param N: number of points
        :param x0: array of coordinates
        :param f0: array of state variables
        :param eps: threshold
        """

        xmin = np.min(x0)
        xmax = np.max(x0)

        rows, cols = u0.shape

        func_array = []

        func_linear_array = []

        for col in range(cols):

            func_array.append(interpolate.make_interp_spline(x0, u0[:,col], order))

            func_linear_array.append(interpolate.make_interp_spline(x0, u0[:,col], 1))

        x_adapt = set(np.linspace(xmin, xmax, 2 ** M))


        NP_adapt = len(x_adapt)

        for j in range(M+1, N):

            xj = np.linspace(xmin, xmax, 2**j)

            xj = np.array(list(set(xj) - x_adapt))

            dj = np.zeros_like(xj)

            NP_available = NPmax - NP_adapt

            if NP_available == 0:
                break

            for col in range(cols):

                func = func_array[col]
                func_linear = func_linear_array[col]

                uj = func(xj)
                uj_linear = func_linear(xj)

                dj += abs(uj-uj_linear)/uref[col]

            x_add = xj[dj > eps]
            d_add = dj[dj > eps]

            ind_d_add = np.argsort(d_add)[::-1]
            x_add_sorted = x_add[ind_d_add]

            x_adapt = set(list(x_adapt) + list(x_add_sorted[0:NP_available]))

            NP_adapt = len(x_adapt)

        x_adapt = np.sort(np.array(list(x_adapt)))

        f_adapt = np.empty((len(x_adapt), cols))
        f_linear = np.empty((len(x_adapt), cols))

        for col in range(cols):

            f_adapt[:,col] = func_array[col](x_adapt)
            f_linear[:,col] = func_linear_array[col](x_adapt)


        return x_adapt, f_adapt.flatten(), f_linear.flatten()
