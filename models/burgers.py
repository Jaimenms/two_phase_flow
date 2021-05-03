from models.model.model import Model
import numpy as np
from methods.fdm.operations.gradient_hrs import GradientHRS
from methods.fdm.operations.second_gradient import SecondGradient
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from models.model.model_domain import ModelDomain
from models.model.model_variable import ModelVariable, ModelRegionEnum
from models.model.model_parameter import ModelParameter
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum


class Burgers(Model):

    jacobian = None
    iter = None

    def __init__(self, x,
                 lb: ModelParameter,
                 ub: ModelParameter,
                 scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N4,
                 scheme_secondorder: SchemeM2FDMEnum = SchemeM2FDMEnum.CENTRAL_N4,
                 flux_delimiter=FluxDelimiterEnum.CUBISTA,
                 ):
        super().__init__()

        self.x = ModelDomain("x", value=x, unit="m", description="x1 coordinate")
        self.register_domain(self.x)

        self.u = ModelVariable("u-velocity", domains=(self.domains['x'],))
        self.register_variables((self.u,))

        self.register_parameters((lb, ub))

        # Operators
        self.grad_x = GradientHRS(self.x.base_value, axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad2_x = SecondGradient(self.x.base_value, axis=0, scheme=scheme_secondorder)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.parse("u-velocity", y)

        dudt = self.parse("u-velocity", yp)

        res_u = dudt + 0.5 * self.grad_x(u**2, u) - 0.001*self.grad2_x(y)

        eqs = []
        if self.parameters['lb'].base_value is None and self.parameters['ub'].base_value is None:
            res_u = self.apply_regions(res_u, regions=(ModelRegionEnum.CLOSED_CLOSED, ))
            eqs.append(res_u)
        elif self.parameters['lb'].base_value is not None:
            res_u = self.apply_regions(res_u, regions=(ModelRegionEnum.OPEN_CLOSED,))
            eqs.append(res_u)
            lb_x = self.apply_regions(u, regions=(ModelRegionEnum.LOWER, )) - self.parameters['lb'].base_value
            #ub_x = self.apply_regions(u, regions=(ModelRegionEnum.UPPER, )) - self.apply_regions(u, regions=(ModelRegionEnum.UPPER_MINUS_ONE,))
            eqs.append(lb_x)
        elif self.parameters['ub'].base_value is not None:
            res_u = self.apply_regions(res_u, regions=(ModelRegionEnum.CLOSED_OPEN,))
            eqs.append(res_u)
            ub_x = self.apply_regions(u, regions=(ModelRegionEnum.UPPER,)) - self.parameters['ub'].base_value
            # ub_x = self.apply_regions(u, regions=(ModelRegionEnum.UPPER, )) - self.apply_regions(u, regions=(ModelRegionEnum.UPPER_MINUS_ONE,))
            eqs.append(ub_x)
        else:
            res_u = self.apply_regions(res_u, regions=(ModelRegionEnum.OPEN_OPEN, ))
            eqs.append(res_u)

            lb_x = self.apply_regions(u, regions=(ModelRegionEnum.LOWER, )) - self.parameters['lb'].base_value
            #lb_x = self.apply_regions(u, regions=(ModelRegionEnum.LOWER, )) - self.apply_regions(u, regions=(ModelRegionEnum.LOWER_PLUS_ONE, ))
            eqs.append(lb_x)

            ub_x = self.apply_regions(u, regions=(ModelRegionEnum.UPPER, )) - self.parameters['ub'].base_value
            #ub_x = self.apply_regions(u, regions=(ModelRegionEnum.UPPER, )) - self.apply_regions(u, regions=(ModelRegionEnum.UPPER_MINUS_ONE,))
            eqs.append(ub_x)

        res = np.concatenate(eqs, axis=None)

        ires = 0

        return res, ires

    def str_equation(self):

        return "dv/dt + v*dv/dx = 0"

    class Parameters:
        y_LB = None
        y_UB = None
        pass



