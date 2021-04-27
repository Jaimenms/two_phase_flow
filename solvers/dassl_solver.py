from .solver import Solver
from models.model import Model
import numpy as np
import dasslc
import os


class DasslSolver(Solver):

    @staticmethod
    def run(
            model: Model,
            t: np.ndarray,
            y0: np.ndarray,
            yp0=None,
            rtol=1e-6,
            atol=1e-8,
            index=None,
            configuration_file=None,
            display=False
    ):

        if configuration_file is None:
            configuration_file = os.path.join(os.path.dirname(__file__),"dassl_default.dat")

        t, y, yp = dasslc.solve(model, t, y0, yp0, model.Parameters, rtol, atol, index, configuration_file, model.jacobian, display)

        return t, y, yp
