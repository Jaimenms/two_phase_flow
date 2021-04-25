from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table import N8TaylorTable


class DownwindSchemeN8FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTable.lambdas()
