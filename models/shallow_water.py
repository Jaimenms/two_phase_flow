import numpy as np

from methods.fdm.operations.gradient import Gradient
from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from models.model.domain import Domain, Domains
from models.model.equation import Equation, Equations, BoundaryConditionEnum, BoundaryCondition
from models.model.model import Model
from models.model.model_plot_mixin import ModelPlotMixin
from models.model.parameter import Parameters
from models.model.variable import RegionEnum, Variables, Variable


class ShallowWater(Model, ModelPlotMixin):

    def __init__(self,
                 x_domain: Domain,
                 scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N8,
                 scheme_hrs: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2,
                 flux_delimiter=FluxDelimiterEnum.CUBISTA2,
                 ):
        super().__init__()

        v = Variable("v", domains=(x_domain,), unit="m/s")
        h = Variable("h", domains=(x_domain,), unit="m")
        self.variables = Variables((v,h))
        self.parameters = Parameters()
        self.domains = Domains((x_domain,))

        # Operators
        self.grad_hrs_x = GradientHRS(self.domains["x"], axis=0, scheme=scheme_hrs, flux_delimiter=flux_delimiter)
        self.grad_x = Gradient(self.domains["x"], axis=0, scheme=scheme)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        #print(t)

        v = self.variables["v"].parse(y)
        dvdt = self.variables["v"].parse(yp)

        h = self.variables["h"].parse(y)
        dhdt = self.variables["h"].parse(yp)

        res_1 = dvdt + 0.5 * self.grad_hrs_x(v**2, v) + 9.81 * self.grad_x(h)
        res_2 = dhdt + h * self.grad_x(v) + v * self.grad_hrs_x(h, h)

        eq1 = Equation(res_1, regions=(RegionEnum.OPEN_OPEN,))
        eq2 = Equation(res_2, regions=(RegionEnum.CLOSED_CLOSED,))

        bc1 = BoundaryCondition(v, 0.0, kind=BoundaryConditionEnum.NEUMANN, regions_1=(RegionEnum.LOWER,), regions_2=(RegionEnum.LOWER_PLUS_ONE,))
        bc2 = BoundaryCondition(v, 0.0, kind=BoundaryConditionEnum.NEUMANN, regions_1=(RegionEnum.UPPER,), regions_2=(RegionEnum.UPPER_MINUS_ONE,))

        eqs = Equations((eq1, eq2, bc1, bc2))

        res = eqs()

        ires = 0

        return res, ires

    def numerical_jacobian_yprime_block(self, t, y, yp, par):
        if self.jac_y_prime is None:
            self.jac_y_prime = super().numerical_jacobian_yprime_block(t, y, yp, par)
        return self.jac_y_prime

    def jacobian(self, t, y, yp, cj, par):
        return self.numerical_jacobian(t, y, yp, cj, par)
