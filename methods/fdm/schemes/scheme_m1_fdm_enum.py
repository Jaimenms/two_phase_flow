from enum import Enum, unique

from methods.fdm.schemes.central.central_scheme_n2_m1_fdm_weights import CentralSchemeN2M1FDMWeights
from methods.fdm.schemes.central.central_scheme_n4_m1_fdm_weights import CentralSchemeN4M1FDMWeights
from methods.fdm.schemes.central.central_scheme_n6_m1_fdm_weights import CentralSchemeN6M1FDMWeights
from methods.fdm.schemes.central.central_scheme_n8_m1_fdm_weights import CentralSchemeN8M1FDMWeights

from methods.fdm.schemes.upwind.upwind_scheme_n1_m1_fdm_weights import UpwindSchemeN1M1FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n2_m1_fdm_weights import UpwindSchemeN2M1FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n4_m1_fdm_weights import UpwindSchemeN4M1FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n6_m1_fdm_weights import UpwindSchemeN6M1FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n8_m1_fdm_weights import UpwindSchemeN8M1FDMWeights

from methods.fdm.schemes.downwind.downwind_scheme_n1_m1_fdm_weights import DownwindSchemeN1M1FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n2_m1_fdm_weights import DownwindSchemeN2M1FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n4_m1_fdm_weights import DownwindSchemeN4M1FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n6_m1_fdm_weights import DownwindSchemeN6M1FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n8_m1_fdm_weights import DownwindSchemeN8M1FDMWeights


@unique
class SchemeM1FDMEnum(Enum):
    CENTRAL_N2 = (2, CentralSchemeN2M1FDMWeights(), 1)
    CENTRAL_N4 = (4, CentralSchemeN4M1FDMWeights(), 1)
    CENTRAL_N6 = (6, CentralSchemeN6M1FDMWeights(), 1)
    CENTRAL_N8 = (8, CentralSchemeN8M1FDMWeights(), 1)
    UPWIND_N1 = (1, UpwindSchemeN1M1FDMWeights(), 1)
    UPWIND_N2 = (2, UpwindSchemeN2M1FDMWeights(), 1)
    UPWIND_N4 = (4, UpwindSchemeN4M1FDMWeights(), 1)
    UPWIND_N6 = (6, UpwindSchemeN6M1FDMWeights(), 1)
    UPWIND_N8 = (8, UpwindSchemeN8M1FDMWeights(), 1)
    DOWNWIND_N1 = (1, DownwindSchemeN1M1FDMWeights(), 1)
    DOWNWIND_N2 = (2, DownwindSchemeN2M1FDMWeights(), 1)
    DOWNWIND_N4 = (4, DownwindSchemeN4M1FDMWeights(), 1)
    DOWNWIND_N6 = (6, DownwindSchemeN6M1FDMWeights(), 1)
    DOWNWIND_N8 = (8, DownwindSchemeN8M1FDMWeights(), 1)