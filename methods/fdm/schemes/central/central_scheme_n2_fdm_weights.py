from methods.fdm.fdm_weights import FDMWeights
from .central_scheme_fdm_mixin import CentralSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table import N2TaylorTable


class CentralSchemeN2FDMWeights(CentralSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTable.lambdas()
