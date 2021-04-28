from models.model.model import Model
import numpy as np
from methods.fdm.fdm_mixin import FDMMixin, FDMEnum, FluxDelimiterEnum

class Burgers(Model, FDMMixin):

    jacobian = None
    iter = None

    def __init__(self, x, order: FDMEnum = FDMEnum.CENTRAL_N8, flux_delimiter: FluxDelimiterEnum = None):
        super().__init__()
        self.x = x

        self.grad_operator = self.Gradient(x, order=order, flux_delimiter=flux_delimiter)

    def __call__(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        if self.iter is None:
            self.iter = 0
        else:
            new_iter = np.count_nonzero(par.solver_t <= t)
            if new_iter > self.iter:
                self.iter = new_iter
                print("{}/{}".format(self.iter, len(par.solver_t)))

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



