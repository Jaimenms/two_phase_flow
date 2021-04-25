from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n4_taylor_table import N4TaylorTable


class CentralSchemeN4FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 4

    lambdas = N4TaylorTable.lambdas()
