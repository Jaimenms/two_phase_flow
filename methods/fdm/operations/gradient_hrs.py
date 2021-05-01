import numpy as np
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.order_fdm_enum import OrderFDMEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum
from methods.fdm.fdm_error import FDMError
from methods.fdm.operations.operation import Operation
from scipy import sparse


class GradientHRS(Operation):

    def __init__(self,
         x,
         order: OrderFDMEnum = OrderFDMEnum.FIRST,
         scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2,
         flux_delimiter: FluxDelimiterEnum = FluxDelimiterEnum.SMART,
         axis=0,
    ):

        if order.value[0] != 1:
            raise FDMError('Flux delimiter is not compatible with this gradient order.')

        x_hrs, id_nodes, id_faces, x_faces, factor_faces = self.create_hrs(x)

        L_hrs = len(x_hrs)
        L = len(x)

        i_p_pos, i_u_pos, i_d_pos = self.calculate_pud_index(L, +1)
        i_p_neg, i_u_neg, i_d_neg = self.calculate_pud_index(L, -1)

        N, scheme_class, M = scheme.value

        if N > L_hrs - 1:
            raise FDMError('Number of nodes is insufficient for the selected gradient order ')

        weights_lrs = scheme_class.matrix(x)
        weights_hrs = scheme_class.matrix(x_hrs)
        weights = np.zeros((L, L_hrs))
        for j, i in enumerate(id_nodes):
            if i in (id_nodes[0], id_nodes[-1]):
                weights[j,0::2] = weights_lrs[j]
            else:
                weights[j, 0:] = weights_hrs[i]

        self.sparse_weights = sparse.csr_matrix(weights)
        self.x = x
        self.x_faces = x_faces
        self.factor_faces = factor_faces
        self.axis = axis
        self.L = L
        self.N = N
        self.x_hrs = x_hrs
        self.L_hrs = L_hrs
        self.id_nodes = id_nodes
        self.flux_delimiter = flux_delimiter.value
        self.i_p_pos=i_p_pos
        self.i_u_pos=i_u_pos
        self.i_d_pos=i_d_pos
        self.i_p_neg=i_p_neg
        self.i_u_neg=i_u_neg
        self.i_d_neg=i_d_neg

    def create_hrs(self, x):

        L = len(x)
        x_hrs = np.empty(2 * L - 1)
        id_nodes = []
        id_faces = []

        id_faces.append(0)
        i = 0
        for i in range(L - 1):
            x_hrs[2 * i] = x[i]
            id_nodes.append(2 * i)
            x_hrs[2 * i + 1] = 0.5 * (x[i] + x[i + 1])
            id_faces.append(2 * i + 1)
        x_hrs[-1] = x[-1]
        id_nodes.append(2 * i + 2)

        x_i = x[0:-1]
        x_ip1 = x[1:]
        x_faces = 0.5 * (x_i + x_ip1)

        factor_faces = (x_faces - x_i) / (x_ip1 - x_i)

        return x_hrs, id_nodes, id_faces, x_faces, factor_faces

    def __call__(self, f, a=None,):
        if a is None:
            a = f
        return self.recursive_slice_call(f, a)

    def slice_call(self, flux_slice: np.ndarray, a_slice: np.ndarray):
        hrs_flux_slice = self.calculate_hrs_flux(flux_slice, a_slice)
        grad_slice = self.sparse_weights.dot(hrs_flux_slice)
        return grad_slice

    def calculate_direction(self, a):
        a_i = a[0:-1]
        a_ip1 = a[1:]
        return np.sign(a_i + self.factor_faces * (a_ip1 - a_i)).astype(int)

    @staticmethod
    def calculate_pud_index(L, direction):

        i_p = np.arange(0, L - 1, 1) + 0.5 * (1 - direction)
        i_p = i_p.astype(int)
        i_u_aux = i_p - direction
        i_d = i_p + direction
        i_u = np.select(
            condlist=[i_u_aux < 0, i_u_aux < L, i_u_aux >= L, ],
            choicelist=[i_p, i_u_aux, i_p],
        )

        return i_p, i_u, i_d

    def get_pointers(self, ind_pos, ind_neg, pos_pointer, neg_pointer):
        pointer = np.empty(self.L - 1, dtype=int)
        pointer[ind_pos] = pos_pointer[ind_pos]
        pointer[ind_neg] = neg_pointer[ind_neg]
        return pointer

    def calculate_hrs_flux(self, flux, a):

        direction = self.calculate_direction(a)
        ind_pos = direction >= 0
        ind_neg = np.invert(ind_pos)

        i_u = self.get_pointers(ind_pos, ind_neg, self.i_u_pos, self.i_u_neg)
        i_p = self.get_pointers(ind_pos, ind_neg, self.i_p_pos, self.i_p_neg)
        i_d = self.get_pointers(ind_pos, ind_neg, self.i_d_pos, self.i_d_neg)

        flux_p = flux[i_p]
        flux_u = flux[i_u]
        flux_d = flux[i_d]

        xp = self.x[i_p]
        xu = self.x[i_u]
        xd = self.x[i_d]

        phi_p = self.normalize(flux_p, flux_u, flux_d)
        tetha_p = self.normalize(xp, xu, xd)
        if direction[0] >= 0:
            tetha_p[0] = tetha_p[1]
            phi_p[0] = phi_p[1]
        if direction[-1] < 0:
            tetha_p[-1] = tetha_p[-2]
            phi_p[-1] = phi_p[-2]

        tetha_f = self.normalize(self.x_faces, xu, xd)
        phi_f = self.flux_delimiter(phi_p, tetha_f, tetha_p)
        flux_f = flux_u + phi_f * (flux_d - flux_u)

        grid_flux = np.empty(2 * self.L - 1, np.float64)
        grid_flux[0::2] = flux
        grid_flux[1::2] = flux_f

        return grid_flux

    @staticmethod
    def normalize(p, u, d):
        eps = 1e-9
        den = d - u
        #den[den == 0] = eps
        den[np.absolute(den) < eps] = eps
        return (p - u) / den
