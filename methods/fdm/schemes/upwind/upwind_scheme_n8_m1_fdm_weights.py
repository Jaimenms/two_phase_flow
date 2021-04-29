from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table_m1 import N8TaylorTableM1


class UpwindSchemeN8M1FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTableM1.lambdas()
