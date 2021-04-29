from models.model.model import Model
import numpy as np
from methods.fdm.fdm_mixin import FDMMixin, SchemeM1FDMEnum, FluxDelimiterEnum
from models.model.model_domain import ModelDomain
from models.model.model_variable import ModelVariable
from models.model.model_parameter import ModelParameter
from models.model.model_constant import ModelConstant


class AdvancedModelMixin:

    DOMAINS = dict()
    CONSTANTS = dict()
    PARAMETERS = dict()
    VARIABLES = dict()
    INDEXES = dict()

    def register_domain(self, input: ModelDomain):
        self.DOMAINS[input.name] = input

    def register_variable(self, input: ModelVariable):
        self.VARIABLES[input.name] = input

        if not self.INDEXES:
            _from = 0
        else:
            _from = max((index.fini for name, index in self.INDEXES.items()))

        n = 1
        for domain in input.domains:
            n *= len(domain)

        _to = _from + n

        self.INDEXES[input.name] = {"ini": _from, "fini": _to}


    def register_parameter(self, input: ModelParameter):
        self.PARAMETERS[input.name] = input

    def register_constant(self, input: ModelConstant):
        self.CONSTANTS[input.name] = input

    def parse_variable_vector(self, y: np.ndarray, reshape=True):
        parsed_y = dict()
        for name, index in self.INDEXES.items():
            values = y[index.ini:index.fini]
            if reshape:
                variable = self.VARIABLES[name]
                values = np.reshape(values, newshape=variable.shape)
            parsed_y[name] = values
        return parsed_y

    def domain(self, name) -> ModelDomain:
        return self.DOMAINS[name]

    def domains(self):
        return self.DOMAINS


class Burgers(Model, FDMMixin, AdvancedModelMixin):

    jacobian = None

    def __init__(self, x1, x2, scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter: FluxDelimiterEnum = None):
        super().__init__()

        # Register all domains
        self.register_domain(ModelDomain("x1", value=x1, unit="m", description="x1 coordinate"))
        self.register_domain(ModelDomain("x2", value=x2, unit="m", description="x2 coordinate"))

        # Register all constants
        # self.register_constant()

        # Register all parameters
        # self.register_parameter()

        # Register all variables
        self.register_variable(ModelVariable("u-velocity", domains=(self.domain('x1'), self.domain('x2'))))
        self.register_variable(ModelVariable("v-velocity", domains=(self.domain('x1'), self.domain('x2'))))

        # Operators
        self.grad_x1 = self.Gradient(self.domain('x1').base_value, axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x2 = self.Gradient(self.domain('x2').base_value, axis=1, scheme=scheme, flux_delimiter=flux_delimiter)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        y = self.parse_variable_vector(y)
        yp = self.parse_variable_vector(yp)

        u = y["u-velocity"]
        dudt = yp["u-velocity"]

        v = y["v-velocity"]
        dvdt = yp["v-velocity"]

        res_u = dudt + u*self.grad_x1(u, a=u) + v*self.grad_x2(u, a=v)

        res_v = dvdt + v*self.grad_x1(v, a=u) + v*self.grad_x2(v, a=v)

        res = np.concatenate((res_u, res_v), axis=None)

        ires = 0

        return res, ires

    def str_equation(self):
        return "dv/dt + v*dv/dx = 0"

    class Parameters:
        y_LB = None
        y_UB = None
        pass



