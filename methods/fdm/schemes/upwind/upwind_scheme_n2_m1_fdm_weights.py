from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table_m1 import N2TaylorTableM1


class UpwindSchemeN2M1FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTableM1.lambdas()
