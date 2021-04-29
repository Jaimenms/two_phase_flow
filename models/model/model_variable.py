import numpy as np
from typing import Tuple, Union
from pint import UnitRegistry
from model_domain import ModelDomain

ureg = UnitRegistry()
Q_ = ureg.Quantity


class ModelVariable:

    def __init__(self, name, domains: Tuple[ModelDomain, ...] = (), description=""):

        self.name = name
        self.domains = domains
        self.description = description
        self.shape = (len(domain) for domain in domains)
