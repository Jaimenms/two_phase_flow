import numpy as np

from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from methods.fdm.operations.second_gradient import SecondGradient, SchemeM2FDMEnum
from models.model.model import Model
from models.model.domain import Domain, Domains
from models.model.equation import Equation, Equations, BoundaryConditionEnum, BoundaryCondition
from models.model.model_plot_mixin import ModelPlotMixin
from models.model.parameter import ConstantParameter, Parameters
from models.model.variable import RegionEnum, Variable, Variables


class Burgers(Model, ModelPlotMixin):

    def __init__(self,
                 x_domain: Domain,
                 scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N4,
                 scheme_secondorder: SchemeM2FDMEnum = SchemeM2FDMEnum.CENTRAL_N4,
                 flux_delimiter=FluxDelimiterEnum.CUBISTA,
                 ):

        u = Variable(name="u", unit="m/s", domains=(x_domain,))
        self.variables = Variables((u,))

        visc = ConstantParameter(name='visc', unit="m**2/s")
        lb = ConstantParameter(name='lb', unit="m/s")
        ub = ConstantParameter(name='ub', unit="m/s")
        self.parameters = Parameters((visc, lb, ub))

        self.domains = Domains((x_domain,))
        self.jac_y_prime = None

        # Operators
        self.grad_x = GradientHRS(self.domains["x"], axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad2_x = SecondGradient(self.domains["x"], axis=0, scheme=scheme_secondorder)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.variables["u"].parse(y)
        dudt = self.variables["u"].parse(yp)

        visc = self.parameters['visc']
        lb = self.parameters['lb']
        ub = self.parameters['ub']

        res_u = dudt + 0.5 * self.grad_x(u**2, u) - visc()*self.grad2_x(y)

        eq_list = []
        if lb() is not None and ub() is not None:
            bc1 = BoundaryCondition(u, lb(), kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,),)
            bc2 = BoundaryCondition(u, ub(), kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER,),)
            eq1 = Equation(res_u, regions=(RegionEnum.OPEN_OPEN,))
            eq_list.append(bc1)
            eq_list.append(bc2)
            eq_list.append(eq1)
        elif lb() is not None:
            bc1 = BoundaryCondition(u, lb(), kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,),)
            eq1 = Equation(res_u, regions=(RegionEnum.OPEN_CLOSED,))
            eq_list.append(bc1)
            eq_list.append(eq1)
        elif ub() is not None:
            bc2 = BoundaryCondition(u, ub(), kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER,),)
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

    def numerical_jacobian_yprime_block(self, t, y, yp, par):
        if self.jac_y_prime is None:
            self.jac_y_prime = super().numerical_jacobian_yprime_block(t, y, yp, par)
        return self.jac_y_prime

    def jacobian(self, t, y, yp, cj, par):
        return self.numerical_jacobian(t, y, yp, cj, par)
