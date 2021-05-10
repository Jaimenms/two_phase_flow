import numpy as np

from methods.fdm.operations.gradient import Gradient
from methods.fdm.operations.gradient_hrs import GradientHRS, SchemeM1FDMEnum, FluxDelimiterEnum
from models.model.domain import Domain
from models.model.equation import Equation, Equations, BoundaryConditionEnum, BoundaryCondition
from models.model.model import Model, Domains, Variables, Parameters
from models.model.model_plot_mixin import ModelPlotMixin
from models.model.parameter import TimeDependentParameter, ConstantParameter
from models.model.variable import RegionEnum
from models.model.variable import Variable
from models.toolbox.dimensionless import Dimensionless
from models.toolbox.geometry import Geometry
from models.toolbox.hydraulics import Hydraulics


class TwoPhaseFlow(Model, ModelPlotMixin):

    def __init__(
            self,
            x_domain: Domain,
            scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N6,
            flux_delimiter=FluxDelimiterEnum.CUBISTA,
    ):

        self.domains = Domains((x_domain,))

        qL = Variable("qL", domains=(x_domain,), unit="kg/s")
        qG = Variable("qG", domains=(x_domain,), unit="kg/s")
        alphaL = Variable("alphaL", domains=(x_domain,), unit="")
        P = Variable("P", domains=(x_domain,), unit="Pa")
        self.variables = Variables((qL, qG, alphaL, P))

        D = ConstantParameter('D', "in")
        epw = ConstantParameter('epw', "um")
        g = ConstantParameter('g', "m/s**2")
        rhoL = ConstantParameter('rhoL', "kg/m**3")
        muL = ConstantParameter('muL', "Pa*s")
        drhoLdP = ConstantParameter('drhoLdP', "kg/m**3/Pa")
        rhoG = ConstantParameter('rhoG', "kg/m**3")
        muG = ConstantParameter('muG', "Pa*s")
        drhoGdP = ConstantParameter('drhoGdP', "kg/m**3/Pa")
        z = ConstantParameter('z', "m")
        qL_lb = TimeDependentParameter('qL_lb', "kg/s")
        qG_lb = TimeDependentParameter('qG_lb', "kg/s")
        P_ub = TimeDependentParameter('P_ub', "Pa")

        self.parameters = Parameters((D, epw, g, rhoL, muL, drhoLdP, rhoG, muG, drhoGdP, z, qL_lb, qG_lb, P_ub))

        # Operators
        self.grad_x_hrs = GradientHRS(self.domains["x"], axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x = Gradient(self.domains["x"], axis=0, scheme=scheme)

    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        D = self.parameters['D']()
        epw = self.parameters['epw']()
        rhoL = self.parameters['rhoL']()
        muL = self.parameters['muL']()
        drhoLdP = self.parameters['drhoLdP']()
        rhoG = self.parameters['rhoG']()
        muG = self.parameters['muG']()
        drhoGdP = self.parameters['drhoGdP']()
        z = self.parameters['z']()
        g = self.parameters['g']()

        qL_lb = self.parameters['qL_lb'](t)
        qG_lb = self.parameters['qG_lb'](t)
        P_ub = self.parameters['P_ub'](t)

        qL = self.variables["qL"].parse(y)
        dqLdt = self.variables["qL"].parse(yp)

        qG = self.variables["qG"].parse(y)
        dqGdt = self.variables["qG"].parse(yp)

        alphaL = self.variables["alphaL"].parse(y)
        dalphaLdt = self.variables["alphaL"].parse(yp)

        P = self.variables["P"].parse(y)
        dPdt = self.variables["P"].parse(yp)

        alphaG = 1-alphaL
        dalphaGdt = -dalphaLdt

        A = Geometry.area(D)
        betha = Geometry.stratified_angle(alphaL)
        perL, perG, perI = Geometry.stratified_perimeters(D, betha)

        vL = Hydraulics.velocity(qL, rhoL, alphaL*A)
        vG = Hydraulics.velocity(qG, rhoG, alphaG*A)

        ReL = Dimensionless.reynolds(D, vL, rhoL, muL)
        ReG = Dimensionless.reynolds(D, vG, rhoG, muG)

        _, fDL = Hydraulics.ff_via_churchil(ReL, epw, D)
        _, fDG = Hydraulics.ff_via_churchil(ReG, epw, D)

        tauL = Hydraulics.shear_stress(fDL, rhoL, vL)
        tauG = Hydraulics.shear_stress(fDG, rhoG, vG)
        tauI = Hydraulics.shear_stress(fDG, rhoG, vG-vL)

        res_1 = rhoL*dalphaLdt + alphaL*drhoLdP*dPdt + 1/A * self.grad_x(qL)
        res_2 = rhoG*dalphaGdt + alphaG*drhoGdP*dPdt + 1/A * self.grad_x(qG)
        res_3 = dqLdt + alphaL * A * self.grad_x(P) + self.grad_x_hrs(qL * vL, qL) + tauL*perL - tauI*perI + alphaL * rhoL * g * A * self.grad_x(z)
        res_4 = dqGdt + alphaG * A * self.grad_x(P) + self.grad_x_hrs(qG * vG, qG) + tauG*perG + tauI*perI + alphaG * rhoG * g * A * self.grad_x(z)

        eq1 = Equation(res_1, regions=(RegionEnum.OPEN_CLOSED,))
        eq2 = Equation(res_2, regions=(RegionEnum.OPEN_CLOSED,))
        eq3 = Equation(res_3, regions=(RegionEnum.CLOSED_OPEN,))
        eq4 = Equation(res_4, regions=(RegionEnum.OPEN_CLOSED,))

        bc1 = BoundaryCondition(qL, qL_lb, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,))
        bc2 = BoundaryCondition(qG, qG_lb, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.LOWER,))
        bc3 = BoundaryCondition(P, P_ub, kind=BoundaryConditionEnum.DIRICHLET, regions_1=(RegionEnum.UPPER,))
        bc4 = BoundaryCondition(alphaL, 0.0, kind=BoundaryConditionEnum.NEUMANN, regions_1=(RegionEnum.LOWER,), regions_2=(RegionEnum.LOWER_PLUS_ONE,))

        eqs = Equations((eq1, eq2, eq3, eq4, bc1, bc2, bc3, bc4))

        res = eqs()

        ires = 0

        return res, ires
