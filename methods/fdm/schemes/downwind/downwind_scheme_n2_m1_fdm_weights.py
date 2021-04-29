from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table_m1 import N2TaylorTableM1


class DownwindSchemeN2M1FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTableM1.lambdas()
