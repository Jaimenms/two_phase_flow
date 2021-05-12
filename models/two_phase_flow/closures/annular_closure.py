import numpy as np

from models.toolbox.dimensionless import Dimensionless
from models.toolbox.geometry import Geometry
from models.toolbox.hydraulics import Hydraulics


class AnnularClosure:

    @staticmethod
    def model_closure(obj, t, qL, qG, alphaL, P):

        D = obj.parameters['D']()
        g = obj.parameters['g']()
        epw = obj.parameters['epw']()
        rhoL = obj.parameters['rhoL']()
        muL = obj.parameters['muL']()
        rhoG = obj.parameters['rhoG']()
        muG = obj.parameters['muG']()

        alphaG = 1-alphaL

        A = Geometry.area(D)

        perG = 0
        perL = np.pi * D
        perI = np.pi * D * alphaG ** 0.5

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

        dPL = 0.0
        dPG = 0.0
        dPI = 0.0

        return tauL*perL, tauG*perG, tauI*perI, dPL, dPG, dPI