from pint import UnitRegistry

ureg = UnitRegistry()


class ModelConstant:

    def __init__(self, name: str, value: float, unit, description=""):
        eng_value = value * ureg(unit)
        self.name = name
        self.eng_value = eng_value
        self.base_value = eng_value.to_base_units().magnitude
        self.base_unit = str(eng_value.to_base_units().units)
        self.description = description

    def __get__(self):
        return self.base_value

    def __str__(self):
        return "<ModelConstant({}, '{}}')>".format(self.base_value, self.base_unit)

