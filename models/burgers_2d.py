import numpy as np

from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from methods.fdm.operations.second_gradient import SecondGradient, SchemeM2FDMEnum
from models.model.model import Model, Domains, Variables, Parameters
from models.model.equation import Equation, Equations
from models.model.variable import RegionEnum
from models.model.boundary_condition import BoundaryCondition, BoundaryConditionEnum
from models.model.model_plot_mixin import ModelPlotMixin


class Burgers2D(Model, ModelPlotMixin):

    def __init__(
            self,
             domains: Domains = Domains(),
             variables: Variables = Variables(),
             parameters: Parameters = Parameters(),
             scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N4,
             scheme_second: SchemeM2FDMEnum = SchemeM2FDMEnum.CENTRAL_N4,
             flux_delimiter: FluxDelimiterEnum = FluxDelimiterEnum.CUBISTA2
        ):
        super().__init__(domains=domains, parameters=parameters, variables=variables)

        self.parameters = parameters
        self.domains = domains
        self.variables = variables

        # Operators
        self.grad_x1 = GradientHRS(self.domains["x1"](), axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x2 = GradientHRS(self.domains["x2"](), axis=1, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad2_x1 = SecondGradient(self.domains["x1"](), axis=0, scheme=scheme_second,)
        self.grad2_x2 = SecondGradient(self.domains["x2"](), axis=1, scheme=scheme_second,)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.variables["u"].parse(y)
        dudt = self.variables["u"].parse(yp)
        visc = self.parameters['visc']()

        res_u = dudt + 0.5*self.grad_x1(u**2, a=u) + 0.5*self.grad_x2(u**2, a=u) + visc*(self.grad2_x1(u) + self.grad2_x2(u))

        eq1 = Equation(res_u, regions=(RegionEnum.OPEN_OPEN,RegionEnum.OPEN_OPEN,))

        bc1 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER, RegionEnum.ALL), )
        bc2 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER, RegionEnum.ALL), )
        bc3 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.OPEN_OPEN, RegionEnum.LOWER), )
        bc4 = BoundaryCondition(u, 0.0, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.OPEN_OPEN, RegionEnum.UPPER), )

        eqs = Equations((eq1, bc1, bc2, bc3, bc4, ))

        res = eqs()

        ires = 0

        return res, ires
