from models.model.model import Model
import numpy as np
from methods.fdm.operations.gradient_hrs import GradientHRS
from methods.fdm.operations.gradient import Gradient
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from models.model.domain import Domains
from models.model.equation import Equation, Equations
from models.model.variable import RegionEnum, Variables
from models.model.parameter import Parameters
from models.model.boundary_condition import BoundaryCondition, BoundaryConditionEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum


class ShallowWater(Model):

    jacobian = None
    iter = None

    def __init__(self,
                 domains: Domains = Domains(),
                 variables: Variables = Variables(),
                 parameters: Parameters = Parameters(),
                 scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N8,
                 scheme_hrs: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2,
                 flux_delimiter=FluxDelimiterEnum.CUBISTA2,
                 ):
        super().__init__()

        self.parameters = parameters
        self.domains = domains
        self.variables = variables

        # Operators
        self.grad_hrs_x = GradientHRS(self.domains["x"](), axis=0, scheme=scheme_hrs, flux_delimiter=flux_delimiter)
        self.grad_x = Gradient(self.domains["x"](), axis=0, scheme=scheme)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        v = self.variables["v"].parse(y)
        dvdt = self.variables["v"].parse(yp)

        h = self.variables["h"].parse(y)
        dhdt = self.variables["h"].parse(yp)

        res_1 = dvdt + 0.5 * self.grad_hrs_x(v**2, v) + 9.81 * self.grad_x(h)
        res_2 = dhdt + self.grad_hrs_x(v*h, v)

        eq1 = Equation(res_1, regions=(RegionEnum.OPEN_OPEN,))
        eq2 = Equation(res_2, regions=(RegionEnum.CLOSED_CLOSED,))

        bc1 = BoundaryCondition(v, 0.0, kind=BoundaryConditionEnum.NEUMANN, regions_1=(RegionEnum.LOWER,), regions_2=(RegionEnum.LOWER_PLUS_ONE,))
        bc2 = BoundaryCondition(v, 0.0, kind=BoundaryConditionEnum.NEUMANN, regions_1=(RegionEnum.UPPER,), regions_2=(RegionEnum.UPPER_MINUS_ONE,))

        eqs = Equations((eq1, eq2, bc1, bc2))

        res = eqs()

        ires = 0

        return res, ires

    def str_equation(self):

        return "..."

    class Parameters:
        pass



