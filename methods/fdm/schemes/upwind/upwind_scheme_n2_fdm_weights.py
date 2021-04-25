from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n2_taylor_table import N2TaylorTable


class UpwindSchemeN2FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N2TaylorTable.lambdas()
