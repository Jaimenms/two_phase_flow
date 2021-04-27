from .model import Model
import numpy as np
from methods.fdm.fdm_mixin import FDMMixin, FDMEnum, FluxDelimiterEnum


class Burgers(Model, FDMMixin):

    jacobian = None

    def __init__(self, x, order: FDMEnum = FDMEnum.CENTRAL_N8, flux_delimiter: FluxDelimiterEnum = None):
        super().__init__()
        self.x = x

        self.grad_operator = self.Gradient(x, order=order, flux_delimiter=flux_delimiter)

    def __call__(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        res = yp + 0.5 * self.grad_operator(y**2, y)
        ires = 0

        return res, ires

    def str_equation(self):

        return "dv/dt + v*dv/dx = 0"

    class Parameters:
        pass



