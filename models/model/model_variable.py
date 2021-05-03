import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from models.model.model_domain import ModelDomain
from enum import Enum, unique

ureg = UnitRegistry()
Q_ = ureg.Quantity


@unique
class ModelRegionEnum(Enum):
    CLOSED_CLOSED = np.s_[0:]
    CLOSED_OPEN = np.s_[0:-1]
    OPEN_OPEN = np.s_[1:-1]
    OPEN_CLOSED = np.s_[1:]
    LOWER = np.s_[0]
    UPPER = np.s_[-1]
    ALL = np.s_[:]
    LOWER_PLUS_ONE = np.s_[1]
    UPPER_MINUS_ONE = np.s_[-2]



class ModelVariable:

    def __init__(self, name, domains: Tuple[ModelDomain, ...] = (), description=""):

        self.name = name
        self.domains = domains
        self.description = description
        self.shape = tuple([len(domain) for domain in domains])
        self.offset = None

        size = 1
        for ele in  self.shape:
            size *= ele
        self.size = size

    def register(self, offset):
        self.offset = offset





