from models.model.model import Model
import numpy as np
from methods.fdm.operations.gradient_hrs import GradientHRS
from methods.fdm.operations.gradient import Gradient
from methods.fdm.operations.second_gradient import SecondGradient
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from models.model.model_domain import ModelDomain
from models.model.model_variable import ModelVariable, ModelRegionEnum
from models.model.model_parameter import ModelParameter
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum


class ShallowWater(Model):

    jacobian = None
    iter = None

    def __init__(self, x,
                 scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N8,
                 scheme_hrs: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2,
                 flux_delimiter=FluxDelimiterEnum.CUBISTA2,
                 ):
        super().__init__()

        self.x = ModelDomain("x", value=x, unit="m", description="x1 coordinate")
        self.register_domain(self.x)

        self.v = ModelVariable("v", domains=(self.domains['x'],))
        self.h = ModelVariable("h", domains=(self.domains['x'],))
        self.register_variables((self.v,self.h,))

        # Operators
        self.grad_hrs_x = GradientHRS(self.x.base_value, axis=0, scheme=scheme_hrs, flux_delimiter=flux_delimiter)
        self.grad_x = Gradient(self.x.base_value, axis=0, scheme=scheme)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        v = self.parse("v", y)
        h = self.parse("h", y)
        dvdt = self.parse("v", yp)
        dhdt = self.parse("h", yp)

        res_v = dvdt + 0.5 * self.grad_hrs_x(v**2, v) + 9.81 * self.grad_x(h)
        res_h = dhdt + self.grad_hrs_x(v*h, v)

        eqs = []
        #res_v = self.apply_regions(res_v, regions=(ModelRegionEnum.CLOSED_CLOSED, ))
        #res_v = self.apply_regions(res_v, regions=(ModelRegionEnum.OPEN_CLOSED, ))
        res_v = self.apply_regions(res_v, regions=(ModelRegionEnum.OPEN_OPEN, ))
        #lb_x1 = self.apply_regions(v, regions=(ModelRegionEnum.LOWER,)) - 0.0
        lb_x1 = self.apply_regions(v, regions=(ModelRegionEnum.LOWER,)) - self.apply_regions(v, regions=(ModelRegionEnum.LOWER_PLUS_ONE,))
        #lb_x2 = self.apply_regions(v, regions=(ModelRegionEnum.UPPER,)) - self.apply_regions(v, regions=(ModelRegionEnum.UPPER_MINUS_ONE,))
        lb_x2 = self.apply_regions(v, regions=(ModelRegionEnum.UPPER,)) - self.apply_regions(v, regions=(ModelRegionEnum.UPPER_MINUS_ONE,))
        eqs.append(res_v)
        eqs.append(lb_x1)
        eqs.append(lb_x2)

        res_h = self.apply_regions(res_h, regions=(ModelRegionEnum.CLOSED_CLOSED,))
        eqs.append(res_h)

        res = np.concatenate(eqs, axis=None)

        ires = 0

        return res, ires

    def str_equation(self):

        return "..."

    class Parameters:
        y_LB = None
        y_UB = None
        pass



