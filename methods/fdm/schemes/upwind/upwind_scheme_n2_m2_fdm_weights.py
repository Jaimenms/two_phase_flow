from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table_m2 import N2TaylorTableM2


class UpwindSchemeN2M2FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTableM2.lambdas()
