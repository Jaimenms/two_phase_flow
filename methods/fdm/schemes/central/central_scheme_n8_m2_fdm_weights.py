from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table_m2 import N8TaylorTableM2


class CentralSchemeN8M2FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTableM2.lambdas()
