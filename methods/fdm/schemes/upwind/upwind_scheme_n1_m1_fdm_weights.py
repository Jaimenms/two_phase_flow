from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n1_taylor_table_m1 import N1TaylorTableM1


class UpwindSchemeN1M1FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N1TaylorTableM1.lambdas()
