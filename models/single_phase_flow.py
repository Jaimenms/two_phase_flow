import numpy as np

from methods.fdm.operations.gradient import Gradient
from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from models.model.domain import Domain
from models.model.equation import Equation, Equations, BoundaryConditionEnum, BoundaryCondition
from models.model.model import Model, Domains, Variables
from models.model.model_plot_mixin import ModelPlotMixin
from models.model.parameter import Parameters, ConstantParameter, TimeDependentParameter
from models.model.variable import RegionEnum, Variable
from models.toolbox.dimensionless import Dimensionless
from models.toolbox.geometry import Geometry
from models.toolbox.hydraulics import Hydraulics


class SinglePhaseFlow(Model, ModelPlotMixin):

    def __init__(
             self,
             x_domain: Domain,
             scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N6,
             flux_delimiter=FluxDelimiterEnum.CUBISTA,
    ):

        self.domains = Domains((x_domain,))

        D = ConstantParameter('D', "in")
        g = ConstantParameter('g', "m/s**2")
        epw = ConstantParameter('epw', "um")
        mu = ConstantParameter('mu', "Pa*s")
        rho = ConstantParameter('rho', "kg/m**3")
        drhodP = ConstantParameter('drhodP', "kg/m**3/Pa")
        z = ConstantParameter('z', "m")
        q_lb = TimeDependentParameter('q_lb', "kg/s")
        P_ub = TimeDependentParameter('P_ub', "Pa")

        self.parameters = Parameters((D, epw, g, rho, mu, drhodP, rho, mu, z, q_lb, P_ub))

        q = Variable("q", domains=(x_domain,), unit="kg/s")
        P = Variable("P", domains=(x_domain,), unit="Pa")
        self.variables = Variables((q, P))

        # Operators
        self.grad_x_hrs = GradientHRS(self.domains["x"], axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x = Gradient(self.domains["x"], axis=0, scheme=scheme)

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

        _, fD = Hydraulics.ff_via_churchill(Re, epw, D)

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
