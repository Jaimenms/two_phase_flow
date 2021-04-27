from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n1_taylor_table import N1TaylorTable


class DownwindSchemeN1FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N1TaylorTable.lambdas()
