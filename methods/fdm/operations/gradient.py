import numpy as np
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from methods.fdm.schemes.order_fdm_enum import OrderFDMEnum
from typing import Union
from methods.fdm.fdm_error import FDMError
from methods.fdm.operations.operation import Operation
from scipy import sparse


class Gradient(Operation):

    def __init__(self,
         x,
         order=OrderFDMEnum.FIRST,
         scheme: Union[SchemeM1FDMEnum,SchemeM2FDMEnum] = SchemeM1FDMEnum.CENTRAL_N2,
         axis=0,
    ):

        L = len(x)
        N, scheme_class, M = scheme.value

        if order.value[0] != scheme.value[2]:
            raise FDMError('Gradient order and finite difference scheme are not consistent.')

        if N > L - 1:
            raise FDMError('Number of nodes is insufficient for the selected gradient order ')

        self.weights = sparse.csr_matrix(scheme_class.matrix(x))
        self.x = x
        self.axis = axis
        self.L = L
        self.N = N

    def __call__(self, f):
        return self.recursive_slice_call(f)

    def slice_call(self, flux_slice: np.ndarray):
        return self.weights.dot(flux_slice)
