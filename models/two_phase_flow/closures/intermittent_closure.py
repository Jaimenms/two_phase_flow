import numpy as np

from models.toolbox.dimensionless import Dimensionless
from models.toolbox.geometry import Geometry
from models.toolbox.hydraulics import Hydraulics
from models.two_phase_flow.closures.stratified_closure import StratifiedClosure
from models.two_phase_flow.closures.bubble_closure import BubbleClosure


class IntermittentClosure:

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
        ReLS = Dimensionless.reynolds(D, vLS, rhoL, muL)
        ReM = Dimensionless.reynolds(D, vM, rhoL, muL)

        ffL, fDL = Hydraulics.ff_via_churchill(ReL, epw, D)
        ffLS, fDLS = Hydraulics.ff_via_churchill(ReLS, epw, D)
        ffM, fDM = Hydraulics.ff_via_churchill(ReM, epw, D)

        # Gas Fraction in Slug
        alphaS = 0.058 * (2 * (0.4 * tension / g / (rhoL - rhoG)) ** 0.5 * (2. * ffM / D * vM ** 3) ** (2/5) * (rhoL / tension) ** (3/5) - 0.725) ** 2
        alphaS = np.fmax(alphaS, np.zeros_like(alphaS))
        alphaS = np.fmin(alphaS, np.ones_like(alphaS) - 0.48)
        Rs = 1 - alphaS

        # Slug length
        if D > 0.038:
            ls = np.exp(-3.287 + 4.589 * (np.log(D / 0.0254)) ** 0.5)
        else:
            ls = 18.1 * D

        # Slug frequency
        Eo = np.sqrt(9.81 * D ** 2. * (rhoL - rhoG) / tension)
        Nf = ((D ** 1.5) * (rhoL * (rhoL - rhoG) * g) ** 0.5) / muL
        n = Nf / (260 + 0.85 * Nf)
        wS = vGS / D * (0.0323 / 2. * vLS / vGS * ((rhoL / (rhoL - rhoG)) ** 0.5) * (fDLS ** -0.5) / ( Eo ** 0.2)) ** n

        # Slug Translational Velocity
        C = 1.2
        gD = (9.81 * D) ** 0.5
        vD = 0.35 * gD * np.sin(tetha) + 0.54 * gD * np.cos(tetha)
        vT = C * vM + vD

        # Slug region lengths
        lu = vT / wS
        ls = np.fmin(ls, lu)
        lf = lu - ls
        lf = np.fmax(lf, ls * 1e-3)
        lu = lf + ls

        # Taylor Bubble holdup
        alphaF = lu / lf * alphaG - ls / lf * alphaS
        alphaF = np.fmax(1e-6, alphaF)
        alphaF = min(1 - 1e-6, alphaF)
        Rf = 1 - alphaF

        # Slug Speeds
        uL = (vLS + vT * (alphaG - alphaS)) / Rs
        ub = (vGS - vT * (alphaG - alphaS)) / alphaS

        # Taylor Bubble Speed
        uf = vT - (vT - uL) * Rs / Rf
        uG = vT - (vT - ub) * alphaS / alphaF

        gammaL1, gammaG1, gammaI1, dPL1, dPG1, dPI1 = StratifiedClosure.model_closure(obj, qL*uL/vL, qG*ub/vG, Rf, P)
        gammaL2, gammaG2, gammaI2, dPL2, dPG2, dPI2 = BubbleClosure.model_closure(obj, qL*uf/vL, qG*uG/vG, Rs, P)

        gammaL = (gammaL1*ls + gammaL2*lf)/lu
        gammaG = (gammaG1*ls + gammaG2*lf)/lu
        gammaI = (gammaI1*ls + gammaI2*lf)/lu
        dPL = (dPL1*ls + dPL2*lf)/lu
        dPG = (dPG1*ls + dPG2*lf)/lu
        dPI = (dPI1*ls + dPI2*lf)/lu

        return gammaL, gammaG, gammaI, dPL, dPG, dPI
