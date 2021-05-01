from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from methods.fdm.schemes.order_fdm_enum import OrderFDMEnum
from methods.fdm.operations.gradient import Gradient


class SecondGradient(Gradient):

    def __init__(self,
         x,
         scheme: SchemeM2FDMEnum = SchemeM2FDMEnum.CENTRAL_N4,
         axis=0,
    ):
        super().__init__(x, scheme=scheme, axis=axis, order=OrderFDMEnum.SECOND)