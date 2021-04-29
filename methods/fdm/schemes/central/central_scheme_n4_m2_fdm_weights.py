from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table_m2 import N4TaylorTableM2


class CentralSchemeN4M2FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTableM2.lambdas()
