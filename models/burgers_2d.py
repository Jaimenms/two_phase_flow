import numpy as np

from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from methods.fdm.operations.second_gradient import SecondGradient, SchemeM2FDMEnum
from models.model.domain import Domain, Domains
from models.model.equation import Equation, Equations, BoundaryConditionEnum, BoundaryCondition
from models.model.model import Model
from models.model.model_plot_mixin import ModelPlotMixin
from models.model.parameter import ConstantParameter, Parameters
from models.model.variable import RegionEnum, Variable, Variables


class Burgers2D(Model, ModelPlotMixin):

    def __init__(
            self,
            x1_domain: Domain,
            x2_domain: Domain,
            scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N4,
            scheme_hrs: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N4,
            scheme_second: SchemeM2FDMEnum = SchemeM2FDMEnum.CENTRAL_N4,
            flux_delimiter: FluxDelimiterEnum = FluxDelimiterEnum.CUBISTA2
        ):
        super().__init__()

        self.domains = Domains((x1_domain, x2_domain))

        u = Variable(name="u", unit="m/s", domains=(x1_domain, x2_domain))
        self.variables = Variables((u,))

        visc = ConstantParameter(name='visc', unit="m**2/s")
        lb = ConstantParameter(name='lb', unit="m/s")
        ub = ConstantParameter(name='ub', unit="m/s")
        self.parameters = Parameters((visc, lb, ub))

        # Operators
        self.grad_x1 = GradientHRS(self.domains["x1"], axis=0, scheme=scheme_hrs, flux_delimiter=flux_delimiter)
        self.grad_x2 = GradientHRS(self.domains["x2"], axis=1, scheme=scheme_hrs, flux_delimiter=flux_delimiter)
        self.grad2_x1 = SecondGradient(self.domains["x1"], axis=0, scheme=scheme_second,)
        self.grad2_x2 = SecondGradient(self.domains["x2"], axis=1, scheme=scheme_second,)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.variables["u"].parse(y)
        dudt = self.variables["u"].parse(yp)
        visc = self.parameters['visc']

        res_u = dudt + 0.5*self.grad_x1(u**2, a=u) + 0.5*self.grad_x2(u**2, a=u) - visc()*(self.grad2_x1(u) + self.grad2_x2(u))

        eq1 = Equation(res_u, regions=(RegionEnum.OPEN_OPEN,RegionEnum.OPEN_OPEN,))

        bc1 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER, RegionEnum.ALL), )
        bc2 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER, RegionEnum.ALL), )
        bc3 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.OPEN_OPEN, RegionEnum.LOWER), )
        bc4 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.OPEN_OPEN, RegionEnum.UPPER), )

        eqs = Equations((eq1, bc1, bc2, bc3, bc4, ))

        res = eqs()

        ires = 0

        return res, ires


    def numerical_jacobian_yprime_block(self, t, y, yp, par):
        if self.jac_y_prime is None:
            self.jac_y_prime = super().numerical_jacobian_yprime_block(t, y, yp, par)
        return self.jac_y_prime

    def jacobian(self, t, y, yp, cj, par):
        return self.numerical_jacobian(t, y, yp, cj, par)
