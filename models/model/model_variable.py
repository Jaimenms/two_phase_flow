import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from model_domain import ModelDomain

ureg = UnitRegistry()
Q_ = ureg.Quantity


class ModelVariable:

    def __init__(self, name, description="", domains=Tuple[ModelDomain]):
        pass