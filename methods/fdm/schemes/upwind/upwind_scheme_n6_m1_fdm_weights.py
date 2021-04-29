from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table_m1 import N6TaylorTableM1


class UpwindSchemeN6M1FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTableM1.lambdas()