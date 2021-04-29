from enum import Enum, unique

from methods.fdm.schemes.central.central_scheme_n2_m2_fdm_weights import CentralSchemeN2M2FDMWeights
from methods.fdm.schemes.central.central_scheme_n4_m2_fdm_weights import CentralSchemeN4M2FDMWeights
from methods.fdm.schemes.central.central_scheme_n6_m2_fdm_weights import CentralSchemeN6M2FDMWeights
from methods.fdm.schemes.central.central_scheme_n8_m2_fdm_weights import CentralSchemeN8M2FDMWeights

from methods.fdm.schemes.upwind.upwind_scheme_n2_m2_fdm_weights import UpwindSchemeN2M2FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n4_m2_fdm_weights import UpwindSchemeN4M2FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n6_m2_fdm_weights import UpwindSchemeN6M2FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n8_m2_fdm_weights import UpwindSchemeN8M2FDMWeights

from methods.fdm.schemes.downwind.downwind_scheme_n2_m2_fdm_weights import DownwindSchemeN2M2FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n4_m2_fdm_weights import DownwindSchemeN4M2FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n6_m2_fdm_weights import DownwindSchemeN6M2FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n8_m2_fdm_weights import DownwindSchemeN8M2FDMWeights


@unique
class SchemeM2FDMEnum(Enum):
    CENTRAL_N4 = (4, CentralSchemeN4M2FDMWeights(), 2)
    CENTRAL_N6 = (6, CentralSchemeN6M2FDMWeights(), 2)
    CENTRAL_N8 = (8, CentralSchemeN8M2FDMWeights(), 2)
    UPWIND_N2 = (2, UpwindSchemeN2M2FDMWeights(), 2)
    UPWIND_N4 = (4, UpwindSchemeN4M2FDMWeights(), 2)
    UPWIND_N6 = (6, UpwindSchemeN6M2FDMWeights(), 2)
    UPWIND_N8 = (8, UpwindSchemeN8M2FDMWeights(), 2)
    DOWNWIND_N2 = (2, DownwindSchemeN2M2FDMWeights(), 2)
    DOWNWIND_N4 = (4, DownwindSchemeN4M2FDMWeights(), 2)
    DOWNWIND_N6 = (6, DownwindSchemeN6M2FDMWeights(), 2)
    DOWNWIND_N8 = (8, DownwindSchemeN8M2FDMWeights(), 2)
