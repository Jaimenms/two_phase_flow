from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table_m1 import N6TaylorTableM1


class CentralSchemeN6M1FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTableM1.lambdas()