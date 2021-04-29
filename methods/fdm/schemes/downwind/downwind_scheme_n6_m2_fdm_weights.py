from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table_m2 import N6TaylorTableM2


class DownwindSchemeN6M2FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTableM2.lambdas()