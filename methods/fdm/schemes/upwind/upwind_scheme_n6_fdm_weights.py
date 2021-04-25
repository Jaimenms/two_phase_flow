from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table import N6TaylorTable


class UpwindSchemeN6FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTable.lambdas()