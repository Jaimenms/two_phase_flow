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

            weights = np.empty((L_grid, N + 1))
            inis = []
            finis = []
            for i in range(L_grid):
                weights[i, :], ini, fini = weight_fcn(L_grid, i, x_grid)
                inis.append(ini)
                finis.append(fini)

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
            j = 0
            for i, (c_i, ini, fini) in enumerate(zip(self.weights, self.inis, self.finis)):
                if i in self.id_nodes:
                    f_i = grid_f[ini:fini]
                    grads[j] = sum(f_i[:] * c_i[:])
                    j += 1

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
            new_x = np.empty(2 * L + 1)
            id_nodes = []
            id_faces = []

            new_x[0] = -0.5*(x[1] - x[0])
            id_faces.append(0)
            for i in range(L - 1):
                new_x[2 * i + 1] = x[i]
                id_nodes.append(2 * i + 1)
                new_x[2 * i + 2] = 0.5 * (x[i] + x[i + 1])
                id_faces.append(2 * i + 2)
            new_x[-2] = x[-1]
            id_nodes.append(2 * i + 3)
            new_x[-1] = x[-1] + 0.5*(x[-1] - x[-2])
            id_faces.append(2 * i + 4)
            return new_x, id_nodes, id_faces

        def merge(self, x, flux, flux_f):
            L = len(x)
            merged_flux = np.empty(self.L_grid)
            for i in range(L):
                merged_flux[2 * i] = flux_f[i]
                merged_flux[2 * i + 1] = flux[i]
            merged_flux[-1] = flux_f[-1]

            return merged_flux

        def get_flux_faces(self, flux, a):
            """

            :param flux:
            :param x:
            :param a: dflux/du
            :return:
            """

            flux_f = np.empty(self.L + 1)
            a_face = np.empty(self.L + 1)
            a_face[0] = a[0]
            a_face[1:-1] = (a[1:] + a[0:-1]) / 2
            a_face[-1] = a[-1]

            for j, a_j in enumerate(a_face):
                if a_j >= 0:
                    i_p = j - 1
                    i_u = i_p - 1
                    i_d = i_p + 1
                else:
                    i_p = j
                    i_u = i_p + 1
                    i_d = i_p - 1

                if i_u < 0 or i_d < 0:
                    flux_f[j] = flux[j]
                elif i_u >= self.L or i_d >= self.L:
                    flux_f[j] = flux[-1]
                else:
                    flux_f[j] = self.calculate_flux_f(flux[i_p], flux[i_u], flux[i_d], self.x[i_p], self.x[i_u], self.x[i_d], self.x[j])

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
