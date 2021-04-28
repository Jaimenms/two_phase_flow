from methods.fdm.fdm_mixin import FDMMixin
from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry

ureg = UnitRegistry()
Q_ = ureg.Quantity
