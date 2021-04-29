from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table_m2 import N4TaylorTableM2


class UpwindSchemeN4M2FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTableM2.lambdas()
