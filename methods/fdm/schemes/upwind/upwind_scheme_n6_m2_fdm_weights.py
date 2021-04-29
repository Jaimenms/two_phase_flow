from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table_m2 import N6TaylorTableM2


class UpwindSchemeN6M2FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTableM2.lambdas()