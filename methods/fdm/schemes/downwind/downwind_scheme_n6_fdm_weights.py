from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table import N6TaylorTable


class DownwindSchemeN6FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTable.lambdas()