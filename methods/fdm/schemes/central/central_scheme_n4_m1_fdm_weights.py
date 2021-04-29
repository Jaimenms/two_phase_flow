from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table_m1 import N4TaylorTableM1


class CentralSchemeN4M1FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTableM1.lambdas()
