import numpy as np
from pint import UnitRegistry
from typing import Tuple, Dict

ureg = UnitRegistry()


class Domain:

    def __init__(self, name, value: np.ndarray = None, unit="", description=""):

        if not ureg(unit).check('[length]'):
            raise TypeError("Domain unit does not have length dimension.")

        eng_value = value * ureg(unit)
        self.name = name
        self.eng_value = eng_value
        self.base_value = eng_value.to_base_units().magnitude
        self.base_unit = str(eng_value.to_base_units().units)
        self.description = description
        self.length=len(value)

    def __call__(self):
        return self.base_value

    def __getitem__(self, index):
        return self.base_value[index]

    def __len__(self):
        return self.length


class Domains:

    def __init__(self, domains: Tuple[Domain, ...] = ()):
        self._offset = 0
        _domains = {}
        for domain in domains:
            _domains[domain.name] = domain

        self.domains = _domains

    def __get__(self) -> Dict[str, Domain]:
        return self.domains

    def __getitem__(self, key) -> Domain:
        return self.domains[key]