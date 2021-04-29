import numpy as np
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from methods.fdm.schemes.order_fdm_enum import OrderFDMEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from scipy import sparse
from typing import Union
from methods.fdm.fdm_error import FDMError
import itertools

class FDMMixin:

    class Gradient:

        def __init__(self,
             x,
             order: OrderFDMEnum = OrderFDMEnum.FIRST,
             scheme: Union[SchemeM1FDMEnum, SchemeM2FDMEnum] = SchemeM1FDMEnum.CENTRAL_N2,
             flux_delimiter: FluxDelimiterEnum = None,
             axis=0,
        ):

            if order.value[0] != scheme.value[2]:
                raise FDMError('Gradient order and finite difference scheme are not consistent.')

            if flux_delimiter is not None and order.value[0] != 1:
                raise FDMError('Flux delimiter is not compatible with this gradient order.')

            if flux_delimiter is None:
                x_grid, id_nodes, id_faces = self.normal_grid(x)
            else:
                x_grid, id_nodes, id_faces = self.hrs_grid(x)

            L_grid = len(x_grid)
            L = len(x)

            N, weight_fcn, M = scheme.value

            if N > L_grid - 1:
                raise FDMError('Number of nodes is insufficient for the selected gradient order ')

            full_weights = np.zeros((L, L_grid))
            for j, i in enumerate(id_nodes):
                if i == id_nodes[0] and flux_delimiter is not None:
                    c_i_ast, ini, fini = weight_fcn(L, j, x)
                    c_i = np.zeros(2*N+1)
                    for k, val in enumerate(c_i_ast):
                        c_i[2*k] = val
                    ini = 0
                    fini = 2*N+1
                elif i == id_nodes[-1] and flux_delimiter is not None:
                    c_i_ast, ini, fini = weight_fcn(L, j, x)
                    c_i = np.zeros(2 * N + 1)
                    for k, val in enumerate(c_i_ast):
                        c_i[2 * k] = val
                    ini = - 2*N - 2
                    fini = -1
                else:
                    c_i, ini, fini = weight_fcn(L_grid, i, x_grid)

                full_weights[j, ini:fini] = c_i

            self.sparse_weights = sparse.csr_matrix(full_weights)

            self.x = x
            self.axis = axis
            self.L = L
            self.N = N
            self.x_grid = x_grid
            self.L_grid = L_grid
            self.id_nodes = id_nodes
            if flux_delimiter is not None:
                self.flux_delimiter = flux_delimiter.value
                self.calculate_grid_flux = self.grid_flux_with_delimiter
            else:
                self.flux_delimiter = None
                self.calculate_grid_flux = self.grid_flux_without_delimiter

        def grid_flux_without_delimiter(self, f, a):
            return f

        def grid_flux_with_delimiter(self, f, a):
            flux_f = self.get_flux_faces(f, a)
            grid_f = self.merge(self.x, f, flux_f)
            return grid_f

        def __call__(self, f, a=None,):

            grads = np.zeros_like(f)

            Ni, Nk = f.shape[:self.axis], f.shape[self.axis + 1:]
            for ii in np.ndindex(Ni):
                for kk in np.ndindex(Nk):
                    grads[ii + np.s_[:, ] + kk] = self.calculate(f[ii + np.s_[:, ] + kk], a[ii + np.s_[:, ] + kk])


            return grads

        def calculate(self, f, a):
            grid_f = self.calculate_grid_flux(f, a)
            grads = self.sparse_weights.dot(grid_f)
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

            x_i = self.x[0:-1]
            x_ip1 = self.x[1:]
            a_i = a[0:-1]
            a_ip1 = a[1:]
            x_f = 0.5 * (x_i + x_ip1)

            a_sign_f = np.sign(a_i + (x_f - x_i) / (x_ip1 - x_i) * (a_ip1 - a_i)).astype(int)

            i_p = np.arange(0, self.L - 1, 1) + 0.5 * ( 1 - a_sign_f)
            i_p = i_p.astype(int)
            i_u_aux = i_p - a_sign_f
            i_d = i_p + a_sign_f
            i_u = np.select(
                condlist=[ i_u_aux < 0, i_u_aux < self.L, i_u_aux >= self.L, ],
                choicelist=[ i_p, i_u_aux, i_p ],
            )

            flux_f = self.calculate_flux_f(flux[i_p], flux[i_u], flux[i_d], self.x[i_p], self.x[i_u], self.x[i_d], x_f, a_sign_f)

            return flux_f

        @staticmethod
        def normalize(p, u, d):
            eps = 1e-9
            den = d - u
            # den[den == 0] = eps
            den[np.absolute(den) < eps] = eps
            return (p - u) / den

        def calculate_flux_f(self, flux_p, flux_u, flux_d, xp, xu, xd, xf, a_f):
            phi_p = self.normalize(flux_p, flux_u, flux_d)
            tetha_p = self.normalize(xp, xu, xd)
            if a_f[0] >= 0:
                tetha_p[0] = tetha_p[1]
            if a_f[-1] < 0:
                tetha_p[-1] = tetha_p[-2]

            tetha_f = self.normalize(xf, xu, xd)
            phi_f = self.flux_delimiter(phi_p, tetha_f, tetha_p)
            return flux_u + phi_f * (flux_d - flux_u)
