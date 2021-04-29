from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table_m2 import N4TaylorTableM2


class DownwindSchemeN4M2FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTableM2.lambdas()
