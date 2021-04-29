from models.model.model import Model
import numpy as np
from methods.fdm.fdm_mixin import FDMMixin, SchemeM1FDMEnum, FluxDelimiterEnum
from models.model.model_domain import ModelDomain
from models.model.model_variable import ModelVariable
from models.model.model_parameter import ModelParameter
from models.model.model_constant import ModelConstant


class AdvancedModelMixin:


    def initialize(self):
        self.domains = dict()
        self.constants = dict()
        self.parameters = dict()
        self.variables = dict()
        self._offset = 0

    def register_domain(self, input: ModelDomain):
        self.domains[input.name] = input

    def register_variable(self, input: ModelVariable):
        input.offset = self._offset
        self._offset += input.size
        self.variables[input.name] = input

    def register_parameter(self, input: ModelParameter):
        self.parameters[input.name] = input

    def register_constant(self, input: ModelConstant):
        self.constants[input.name] = input

    def parse(self, name, y: np.ndarray):
        variable = self.variables[name]
        ini = variable.offset
        fini = variable.offset + variable.size
        return np.reshape(y[ini:fini], newshape=variable.shape)

class Burgers(Model, FDMMixin, AdvancedModelMixin):

    jacobian = None

    def __init__(self, x1, x2, scheme: SchemeM1FDMEnum = SchemeM1FDMEnum.CENTRAL_N2, flux_delimiter: FluxDelimiterEnum = None):
        super().__init__()
        self.initialize()

        # Register all domains
        self.register_domain(ModelDomain("x1", value=x1, unit="m", description="x1 coordinate"))
        self.register_domain(ModelDomain("x2", value=x2, unit="m", description="x2 coordinate"))

        # Register all constants
        # self.register_constant()

        # Register all parameters
        # self.register_parameter()

        # Register all variables
        self.register_variable(ModelVariable("u-velocity", domains=(self.domains['x1'], self.domains['x2'])))
        self.register_variable(ModelVariable("v-velocity", domains=(self.domains['x1'], self.domains['x2'])))

        # Operators
        self.grad_x1 = self.Gradient(self.domains['x1'].base_value, axis=0, scheme=scheme, flux_delimiter=flux_delimiter)
        self.grad_x2 = self.Gradient(self.domains['x2'].base_value, axis=1, scheme=scheme, flux_delimiter=flux_delimiter)


    def residue(self, t: float, y: np.ndarray, yp: np.ndarray, par=None):

        u = self.parse("u-velocity", y)
        v = self.parse("v-velocity", y)

        dudt = self.parse("u-velocity", yp)
        dvdt = self.parse("v-velocity", yp)

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



