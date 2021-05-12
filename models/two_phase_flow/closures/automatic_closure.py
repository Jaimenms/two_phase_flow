import numpy as np

from models.toolbox.dimensionless import Dimensionless
from models.toolbox.geometry import Geometry
from models.toolbox.hydraulics import Hydraulics


class AutomaticClosure:

    @staticmethod
    def model_closure(obj, t, qL, qG, alphaL, P):

        D = obj.parameters['D']()
        g = obj.parameters['g']()
        epw = obj.parameters['epw']()
        rhoL = obj.parameters['rhoL']()
        muL = obj.parameters['muL']()
        rhoG = obj.parameters['rhoG']()
        muG = obj.parameters['muG']()
        tetha = obj.parameters['tetha']()

        alphaG = 1-alphaL

        A = Geometry.area(D)
        betha = Geometry.stratified_angle(alphaL)
        hD = 0.5 * (1 - np.cos(0.5 * betha))
        perL, perG, perI = Geometry.stratified_perimeters(D, betha)

        vL = Hydraulics.velocity(qL, rhoL, alphaL*A)
        vG = Hydraulics.velocity(qG, rhoG, alphaG*A)

        DL = 4. * alphaL * A / perL
        DG = 4. * alphaG * A / (perG + perI)

        ReL = Dimensionless.reynolds(DL, vL, rhoL, muL)
        ReG = Dimensionless.reynolds(DG, vG, rhoG, muG)

        _, fDL = Hydraulics.ff_via_churchill(ReL, epw, DL)
        _, fDG = Hydraulics.ff_via_churchill(ReG, epw, DG)

        tauL = Hydraulics.shear_stress(fDL, rhoL, vL)
        tauG = Hydraulics.shear_stress(fDG, rhoG, vG)
        tauI = Hydraulics.shear_stress(fDG, rhoG, vG-vL)

        zcml = D * (0.5 - 1 / 3 / np.pi / alphaL * np.sin(betha / 2) ** 3)
        zcmg = D * (0.5 + 1 / 3 / np.pi / alphaG * np.sin(betha / 2) ** 3)

        dPLG = ((rhoL - rhoG) * hD * D - rhoL * zcml + rhoG * zcmg) * g * np.cos(tetha)
        dPL = alphaG * dPLG
        dPG = -alphaL * dPLG
        dPI = (alphaL * rhoL * (zcml - hD * D) + alphaG * rhoG * (zcmg - hD * D)) * g * np.cos(tetha)

        return tauL*perL, tauG*perG, tauI*perI, dPL, dPG, dPI