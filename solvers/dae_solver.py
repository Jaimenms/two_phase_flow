from .solver import Solver
from models.model import Model
import numpy as np
import dasslc


class DAESolver(Solver):

    def run(self, model: Model, t: np.ndarray, y0: np.ndarray, yp0=None, par=None, rtol=1e-6, atol=1e-8, index=None):

        t, y, yp = dasslc.solve(model, t, y0, yp0, par, rtol, atol, index)

        return t, y, yp
