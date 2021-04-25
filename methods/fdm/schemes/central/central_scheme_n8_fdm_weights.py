from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n8_taylor_table import N8TaylorTable


class CentralSchemeN8FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 8

    lambdas = N8TaylorTable.lambdas()
