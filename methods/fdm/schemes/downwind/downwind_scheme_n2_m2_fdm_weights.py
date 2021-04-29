from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table_m2 import N2TaylorTableM2


class DownwindSchemeN2M2FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTableM2.lambdas()
