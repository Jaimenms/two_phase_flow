from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table_m2 import N8TaylorTableM2


class UpwindSchemeN8M2FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTableM2.lambdas()
