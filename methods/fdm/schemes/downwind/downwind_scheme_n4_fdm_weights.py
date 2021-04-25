from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table import N4TaylorTable


class DownwindSchemeN4FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTable.lambdas()
