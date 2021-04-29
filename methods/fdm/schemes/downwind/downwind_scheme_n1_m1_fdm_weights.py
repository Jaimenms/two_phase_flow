from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n1_taylor_table_m1 import N1TaylorTableM1


class DownwindSchemeN1M1FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N1TaylorTableM1.lambdas()
