from enum import unique, Enum

import numpy as np
from typing import Tuple, Union
from models.model.variable import RegionEnum

@unique
class BoundaryConditionEnum(Enum):
    DIRICHLET = 1
    NEUMANN = 2


class BoundaryCondition:

    def __init__(self,
                 variable_array: np.ndarray,
                 bc_value: Union[float,np.ndarray],
                 kind: BoundaryConditionEnum,
                 regions_1: Tuple[RegionEnum,...],
                 regions_2: Tuple[RegionEnum,...]=(),
                 ):
        if kind == BoundaryConditionEnum.DIRICHLET:
            slices_1 = tuple([region.value for region in regions_1])
            val_1 = variable_array[slices_1]
            self.res = val_1 - bc_value
        else:
            slices_1 = tuple([region.value for region in regions_1])
            slices_2 = tuple([region.value for region in regions_2])
            val_1 = variable_array[slices_1]
            val_2 = variable_array[slices_2]
            self.res = val_1 - val_2 - bc_value

    def __call__(self, *args, **kwargs):
        return self.res


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
