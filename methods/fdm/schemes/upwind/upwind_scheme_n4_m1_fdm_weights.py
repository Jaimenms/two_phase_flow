from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table_m1 import N4TaylorTableM1


class UpwindSchemeN4M1FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTableM1.lambdas()
