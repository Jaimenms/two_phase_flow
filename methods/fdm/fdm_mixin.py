import numpy as np
from methods.fdm.fdm_enum import FDMEnum
from methods.fdm.fdm_error import FDMError
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum


class FDMMixin:

    class Gradient:

        def __init__(self,
                     x: np.ndarray = None,
                     order: FDMEnum = FDMEnum.CENTRAL_N2,
                     flux_delimiter: FluxDelimiterEnum = None
                     ):

            if flux_delimiter is None:
                x_grid, id_nodes, id_faces = self.normal_grid(x)
            else:
                x_grid, id_nodes, id_faces = self.hrs_grid(x)

            L_grid = len(x_grid)
            L = len(x)

            N, weight_fcn = order.value

            if N > L_grid - 1:
                raise FDMError('Number of nodes is insufficient for the selected gradient order ')

            weights = np.empty((L, N + 1))
            inis = []
            finis = []
            j=0
            for i in range(L_grid):
                if i in id_nodes:
                    weights[j, :], ini, fini = weight_fcn(L_grid, i, x_grid)
                    inis.append(ini)
                    finis.append(fini)
                    j += 1

            self.x = x
            self.L = L
            self.N = N
            self.x_grid = x_grid
            self.L_grid = L_grid
            self.weights = weights
            self.inis = inis
            self.finis = finis
            self.id_nodes = id_nodes
            if flux_delimiter is not None:
                self.flux_delimiter = flux_delimiter.value
            else:
                self.flux_delimiter = None

        def __call__(self, f, a=None):

            if a is None:
                a = np.ones_like(f)

            if self.flux_delimiter is not None:
                flux_f = self.get_flux_faces(f, a)
                grid_f = self.merge(self.x, f, flux_f)
            else:
                grid_f = f

            grads = np.empty_like(f)
            for i, (c_i, ini, fini) in enumerate(zip(self.weights, self.inis, self.finis)):
                f_i = grid_f[ini:fini]

                if self.flux_delimiter is not None:
                    if i == 0:
                        c_i = np.zeros_like(c_i)
                        c_i[0] = -1/(self.x[2]-self.x[0])
                        c_i[2] = 1/(self.x[2] - self.x[0])
                    elif i == (self.L-1):
                        c_i = np.zeros_like(c_i)
                        c_i[0] = -1/(self.x[-1]-self.x[-3])
                        c_i[2] = 1/(self.x[-1] - self.x[-3])
                grads[i] = sum(f_i[:] * c_i[:])

            return grads

        def normal_grid(self, x):
            L = len(x)
            new_x = x
            id_nodes = []
            id_faces = []
            for i in range(L):
                id_nodes.append(i)
            return new_x, id_nodes, id_faces

        def hrs_grid(self, x):
            L = len(x)
            new_x = np.empty(2 * L - 1)
            id_nodes = []
            id_faces = []

            id_faces.append(0)
            for i in range(L - 1):
                new_x[2 * i] = x[i]
                id_nodes.append(2 * i)
                new_x[2 * i + 1] = 0.5 * (x[i] + x[i + 1])
                id_faces.append(2 * i + 1)
            new_x[-1] = x[-1]
            id_nodes.append(2 * i + 2)
            return new_x, id_nodes, id_faces

        def merge(self, x, flux, flux_f):
            L = len(x)
            merged_flux = np.empty(self.L_grid)
            for i in range(L - 1):
                merged_flux[2 * i] = flux[i]
                merged_flux[2 * i + 1] = flux_f[i]
            merged_flux[-1] = flux[-1]

            return merged_flux

        def get_flux_faces(self, flux, a):
            """

            :param flux:
            :param x:
            :param a: dflux/du
            :return:
            """

            flux_f = np.empty(self.L - 1)

            for f in range(self.L-1):

                x_f = self.x_grid[2 * f + 1]
                x_i = self.x_grid[2 * f]
                x_ip1 = self.x_grid[2 * f + 2]
                a_i = a[f]
                a_ip1 = a[f + 1]

                a_f = a_i + (x_f - x_i) / (x_ip1 - x_i) * (a_ip1 + a_i)

                if a_f >= 0:
                    i_p = f
                    i_u = i_p - 1
                    i_d = i_p + 1

                else:
                    i_p = f + 1
                    i_u = i_p + 1
                    i_d = i_p - 1

                if i_u < 0 or i_u >= self.L:
                    i_u = i_p

                flux_f[f] = self.calculate_flux_f(flux[i_p], flux[i_u], flux[i_d], self.x[i_p], self.x[i_u], self.x[i_d], self.x[f])

            return flux_f

        def calculate_phi_p(self, flux_p, flux_u, flux_d):
            return (flux_p - flux_u) / (flux_d - flux_u)

        def calculate_tetha_p(self, xp, xu, xd):
            return (xp - xu) / (xd - xu)

        def calculate_tetha_f(self, xf, xu, xd):
            return (xf - xu) / (xd - xu)

        def calculate_flux_f(self, flux_p, flux_u, flux_d, xp, xu, xd, xf):
            if (flux_p - flux_u) == (flux_d - flux_u):
                return flux_d
            phi_p = self.calculate_phi_p(flux_p, flux_u, flux_d)
            tetha_p = self.calculate_tetha_p(xp, xu, xd)
            tetha_f = self.calculate_tetha_f(xf, xu, xd)
            phi_f = self.flux_delimiter(phi_p, tetha_f, tetha_p)
            return flux_u + phi_f * (flux_d - flux_u)
