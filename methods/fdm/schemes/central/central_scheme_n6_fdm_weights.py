from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n6_taylor_table import N6TaylorTable


class CentralSchemeN6FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 6

    lambdas = N6TaylorTable.lambdas()