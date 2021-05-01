from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.order_fdm_enum import OrderFDMEnum
from methods.fdm.operations.gradient import Gradient


class FirstGradient(Gradient):

    def __init__(self,
         x,
         scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2,
         axis=0,
    ):
        super().__init__(x, scheme=scheme, axis=axis, order=OrderFDMEnum.FIRST)