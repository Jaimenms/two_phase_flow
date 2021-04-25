from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table import N2TaylorTable


class DownwindSchemeN2FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTable.lambdas()
