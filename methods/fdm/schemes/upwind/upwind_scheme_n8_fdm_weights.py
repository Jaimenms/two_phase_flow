from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table import N8TaylorTable


class UpwindSchemeN8FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTable.lambdas()
