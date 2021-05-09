from .solver import Solver
from typing import Union
from models.model.model import Model
import numpy as np
import dasslc
import os
from scipy.optimize import root


class DasslSolver(Solver):

    @staticmethod
    def run(
            model: Model,
            t: Union[None, np.ndarray],
            y0: np.ndarray,
            yp0=None,
            rtol=1e-6,
            atol=1e-8,
            index=None,
            configuration_file=None,
            display=False,
            verbose=False
    ):

        if configuration_file is None:
            if verbose:
                configuration_file = os.path.join(os.path.dirname(__file__),"dassl_verbose.dat")
            else:
                configuration_file = os.path.join(os.path.dirname(__file__),"dassl_default.dat")

        if t is None:
            out = root(model.ss_residue, y0, tol=atol)
            y = out['x']
            yp = np.zeros_like(y)
        else:
            t, y, yp = dasslc.solve(model, t, y0, yp0, model.Parameters, rtol, atol, index, configuration_file, model.jacobian, display)

        return t, y, yp
