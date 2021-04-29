from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table_m1 import N2TaylorTableM1


class CentralSchemeN2M1FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTableM1.lambdas()
