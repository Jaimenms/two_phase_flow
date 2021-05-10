from models.model.model import Model
import numpy as np
from methods.fdm.operations.gradient_hrs import GradientHRS
from methods.fdm.operations.second_gradient import SecondGradient
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum
from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from models.model.domain import Domains
from models.model.equation import Equation, Equations
from models.model.variable import RegionEnum, Variables
from models.model.parameter import Parameters
from models.model.boundary_condition import BoundaryCondition, BoundaryConditionEnum
from methods.fdm.flux_delimiters.flux_delimiter_enum import FluxDelimiterEnum


class Burgers(Model):

    jacobian = None
    iter = None

    def __init__(self,
                 domains: Domains = Domains(),
                 variables: Variables = Variables(),
                 parameters: Parameters = Parameters(),
                 scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N4,
                 scheme_secondorder: SchemeM2FDMEnum = SchemeM2FDMEnum.CENTRAL_N4,
                 flux_delimiter=FluxDelimiterEnum.CUBISTA,
                 ):
        super().__init__()

        self.parameters = parameters
        self.domains = domains
        self.variables = variables

        # Operators
        self.grad_x = GradientHRS(self.domains["x"](), axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad2_x = SecondGradient(self.domains["x"](), axis=0, scheme=scheme_secondorder)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.variables["u"].parse(y)
        dudt = self.variables["u"].parse(yp)

        res_u = dudt + 0.5 * self.grad_x(u**2, u) - 0.001*self.grad2_x(y)

        lb = self.parameters['lb']()
        ub = self.parameters['ub']()

        eq_list = []
        if lb is not None and ub is not None:
            bc1 = BoundaryCondition(u, lb, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,),)
            bc2 = BoundaryCondition(u, ub, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER,),)
            eq1 = Equation(res_u, regions=(RegionEnum.OPEN_OPEN,))
            eq_list.append(bc1)
            eq_list.append(bc2)
            eq_list.append(eq1)
        elif lb is not None:
            bc1 = BoundaryCondition(u, lb, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,),)
            eq1 = Equation(res_u, regions=(RegionEnum.OPEN_CLOSED,))
            eq_list.append(bc1)
            eq_list.append(eq1)
        elif ub is not None:
            bc2 = BoundaryCondition(u, ub, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER,),)
            eq1 = Equation(res_u, regions=(RegionEnum.CLOSED_OPEN,))
            eq_list.append(bc2)
            eq_list.append(eq1)
        else:
            eq1 = Equation(res_u, regions=(RegionEnum.CLOSED_CLOSED,))
            eq_list.append(eq1)

        eqs = Equations(tuple(eq_list))

        res = eqs()

        ires = 0

        return res, ires

    def str_equation(self):

        return "dv/dt + v*dv/dx = 0"



