import numpy as np
from typing import Tuple, Union
from models.model.variable import RegionEnum
from models.model.boundary_condition import BoundaryCondition


class Equation:

    def __init__(self, res: np.ndarray, regions: Tuple[RegionEnum,...]):
        slices = tuple([region.value for region in regions])
        self.res = res[slices]

    def __call__(self):
        return self.res



class Equations:

    def __init__(self, equations: Tuple[Union[Equation,BoundaryCondition],...] = ()):
        self.equations = equations

    def __call__(self):
        return np.concatenate([equation() for equation in self.equations], axis=None)