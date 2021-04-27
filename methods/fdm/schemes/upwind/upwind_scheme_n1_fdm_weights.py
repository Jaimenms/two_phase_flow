from methods.fdm.fdm_weights import FDMWeights
from .upwind_scheme_fdm_mixin import UpwindSchemeFDMMixin
from methods.fdm.taylor_tables.n1_taylor_table import N1TaylorTable


class UpwindSchemeN1FDMWeights(UpwindSchemeFDMMixin, FDMWeights):

    N = 2

    lambdas = N1TaylorTable.lambdas()
