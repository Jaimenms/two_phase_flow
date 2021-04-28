from typing import Tuple
from pint import UnitRegistry
from model_domain import ModelDomain

ureg = UnitRegistry()


class ModelParameter:

    def __init__(self, name: str, value: float, unit: str, domains=Tuple[ModelDomain], description=""):

        eng_value = value * ureg(unit)
        self.name = name
        self.eng_value = eng_value
        self.base_value = eng_value.to_base_units().magnitude
        self.base_unit = str(eng_value.to_base_units().units)
        self.domains = domains
        self.description = description

    def __get__(self):
        return self.base_value

    def __str__(self):
        return "<ModelParameter({}, '{}}')>".format(self.base_value, self.base_unit)
