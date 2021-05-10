import numpy as np
from methods.fdm.operations.gradient import Gradient
from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from models.model.model import Model, Domains, Variables, Parameters
from models.model.model_plot_mixin import ModelPlotMixin
from models.model.equation import Equation, Equations
from models.model.variable import RegionEnum
from models.model.boundary_condition import BoundaryCondition, BoundaryConditionEnum
from models.toolbox.dimensionless import Dimensionless
from models.toolbox.hydraulics import Hydraulics
from models.toolbox.geometry import Geometry


class SinglePhaseFlow(Model, ModelPlotMixin):

    def __init__(
             self,
             domains: Domains = Domains(),
             variables: Variables = Variables(),
             parameters: Parameters = Parameters(),
             scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N6,
             flux_delimiter=FluxDelimiterEnum.CUBISTA,
    ):
        super().__init__(domains=domains, parameters=parameters, variables=variables)

        self.parameters = parameters
        self.domains = domains
        self.variables = variables

        # Operators
        self.grad_x_hrs = GradientHRS(self.domains["x"](), axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x = Gradient(self.domains["x"](), axis=0, scheme=scheme)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        D = self.parameters['D']()
        epw = self.parameters['epw']()
        rho = self.parameters['rho']()
        mu = self.parameters['mu']()
        drhodP = self.parameters['drhodP']()
        z = self.parameters['z']()
        g = self.parameters['g']()
        q_lb = self.parameters['q_lb'](t)
        P_ub = self.parameters['P_ub'](t)

        q = self.variables["q"].parse(y)
        dqdt = self.variables["q"].parse(yp)

        P = self.variables["P"].parse(y)
        dPdt = self.variables["P"].parse(yp)

        A = Geometry.area(D)
        per = np.pi * D

        v = Hydraulics.velocity(q, rho, A)

        Re = Dimensionless.reynolds(D, v, rho, mu)

        _, fD = Hydraulics.ff_via_churchil(Re, epw, D)

        tau = Hydraulics.shear_stress(fD, rho, v)

        res_1 = drhodP*dPdt + 1/A * self.grad_x(q)
        res_2 = dqdt + A * self.grad_x(P) + self.grad_x_hrs(q * v, q) + tau * per + rho * g * A * self.grad_x(z)

        eq1 = Equation(res_1, regions=(RegionEnum.OPEN_CLOSED,))
        eq2 = Equation(res_2, regions=(RegionEnum.CLOSED_OPEN,))

        bc1 = BoundaryCondition(q, q_lb, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,))
        bc2 = BoundaryCondition(P, P_ub, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER,))

        eqs = Equations((eq1, eq2, bc1, bc2))

        res = eqs()

        ires = 0

        return res, ires
