from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table import N4TaylorTable


class UpwindSchemeN4FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTable.lambdas()
