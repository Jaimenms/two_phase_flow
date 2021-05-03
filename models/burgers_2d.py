from models.model.model import Model
import numpy as np
from methods.fdm.operations.gradient_hrs import GradientHRS
from methods.fdm.operations.second_gradient import SecondGradient
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from models.model.model_domain import ModelDomain
from models.model.model_variable import ModelVariable, ModelRegionEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum


class Burgers2D(Model,):

    jacobian = None

    def __init__(self, x1, x2, scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter: FluxDelimiterEnum = FluxDelimiterEnum.CUBISTA2):
        super().__init__()

        self.x1 = ModelDomain("x1", value=x1, unit="m", description="x1 coordinate")
        self.x2 = ModelDomain("x2", value=x2, unit="m", description="x2 coordinate")

        # Register all domains
        self.register_domain(self.x1)
        self.register_domain(self.x2)

        # Register all constants
        # self.register_constant()

        # Register all parameters
        # self.register_parameter()

        # Register all variables
        self.u = ModelVariable("u-velocity", domains=(self.domains['x1'], self.domains['x2']))
        self.register_variable(self.u)

        # Operators
        self.grad_x1 = GradientHRS(self.x1.base_value, axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x2 = GradientHRS(self.x2.base_value, axis=1, scheme=scheme, flux_delimiter=flux_delimiter)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.parse("u-velocity", y)

        dudt = self.parse("u-velocity", yp)

        res_u = dudt + 0.5*self.grad_x1(u**2, a=u) + 0.5*self.grad_x2(u**2, a=u)

        res_u = self.apply_regions(res_u, regions=(ModelRegionEnum.OPEN_OPEN, ModelRegionEnum.OPEN_OPEN))

        lb_x1 = self.apply_regions(u, regions=(ModelRegionEnum.LOWER, ModelRegionEnum.ALL)) - 0.0
        ub_x1 = self.apply_regions(u, regions=(ModelRegionEnum.UPPER, ModelRegionEnum.ALL)) - 0.0
        lb_x2 = self.apply_regions(u, regions=(ModelRegionEnum.OPEN_OPEN, ModelRegionEnum.LOWER)) - 0.0
        ub_x2 = self.apply_regions(u, regions=(ModelRegionEnum.OPEN_OPEN, ModelRegionEnum.UPPER)) - 0.0

        res = np.concatenate((res_u, lb_x1, ub_x1, lb_x2, ub_x2), axis=None)

        ires = 0

        return res, ires

    def str_equation(self):
        return "not defined"

    class Parameters:
        y_LB = None
        y_UB = None
        pass



