from methods.fdm.fdm_weights import FDMWeights
from .downwind_scheme_fdm_mixin import DownwindSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table_m2 import N8TaylorTableM2


class DownwindSchemeN8M2FDMWeights(DownwindSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTableM2.lambdas()
