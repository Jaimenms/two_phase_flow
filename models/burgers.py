from models.model.model import Model
import numpy as np
from methods.fdm.operations.gradient_hrs import GradientHRS
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum

class Burgers(Model):

    jacobian = None
    iter = None

    def __init__(self, x, scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter=FluxDelimiterEnum.CUBISTA):
        super().__init__()
        self.x = x

        self.grad_operator = GradientHRS(x, scheme=scheme, flux_delimiter=flux_delimiter)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        res = yp + 0.5 * self.grad_operator(y**2, y)
        if par is not None and par.y_LB is not None:
            res[0] = y[0] - par.y_LB
        if par is not None and par.y_UB is not None:
            res[-1] = y[-1] - par.y_UB
        ires = 0

        return res, ires

    def str_equation(self):

        return "dv/dt + v*dv/dx = 0"

    class Parameters:
        y_LB = None
        y_UB = None
        pass



