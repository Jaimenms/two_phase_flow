from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table_m1 import N8TaylorTableM1


class CentralSchemeN8M1FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTableM1.lambdas()
