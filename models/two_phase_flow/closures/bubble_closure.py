import numpy as np

from models.toolbox.dimensionless import Dimensionless
from models.toolbox.geometry import Geometry
from models.toolbox.hydraulics import Hydraulics


class BubbleClosure:

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
        tension = obj.parameters['tension']()

        alphaG = 1-alphaL

        A = Geometry.area(D)

        vL = Hydraulics.velocity(qL, rhoL, alphaL*A)
        vG = Hydraulics.velocity(qG, rhoG, alphaG*A)

        vLS = Hydraulics.velocity(qL, rhoL, A)
        vGS = Hydraulics.velocity(qG, rhoG, A)

        vM = vGS + vLS

        ReL = Dimensionless.reynolds(D, vL, rhoL, muL)
        ReM = Dimensionless.reynolds(D, vM, rhoL, muL)

        _, fDL = Hydraulics.ff_via_churchill(ReL, epw, D)
        _, fDM = Hydraulics.ff_via_churchill(ReM, epw, D)
        fDG = 0.0

        dbmax = (0.725 + 4.15 * (vGS / vM) ** 0.5) * (tension / rhoL) ** (3 / 5) * ( 2 * fDM * vM ** 3 / D) ** (-2/5)
        db = (16. * tension / g / (rhoL - rhoG)) ** 0.5
        db = np.fmax(db, dbmax)
        vb = (0.51 * db + 2.14 / db) ** 0.5
        CDb = 4 / 3 * ((rhoL - rhoG) / rhoL) * g * db / (vb ** 2)
        fDI = 1.5 * alphaG * rhoL * CDb / db

        perL = np.pi * D * alphaL
        perG = np.pi * D * alphaG
        perI = np.pi * D ** 2 * alphaG / db

        tauL = Hydraulics.shear_stress(fDL, rhoL, vL)
        tauG = Hydraulics.shear_stress(fDG, rhoG, vG)
        tauI = Hydraulics.shear_stress(fDI, rhoG, vG-vL)

        dPL = 1 / 4. * rhoL * (vG - vL) ** 2
        dPG = 0
        dPI = 0

        return tauL*perL, tauG*perG, tauI*perI, dPL, dPG, dPI