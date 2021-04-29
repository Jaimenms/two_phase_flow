import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from models.model.model_domain import ModelDomain

ureg = UnitRegistry()
Q_ = ureg.Quantity


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


