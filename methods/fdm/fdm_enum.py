from enum import Enum, unique

from methods.fdm.schemes.central.central_scheme_n2_fdm_weights import CentralSchemeN2FDMWeights
from methods.fdm.schemes.central.central_scheme_n4_fdm_weights import CentralSchemeN4FDMWeights
from methods.fdm.schemes.central.central_scheme_n6_fdm_weights import CentralSchemeN6FDMWeights
from methods.fdm.schemes.central.central_scheme_n8_fdm_weights import CentralSchemeN8FDMWeights

from methods.fdm.schemes.upwind.upwind_scheme_n2_fdm_weights import UpwindSchemeN2FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n4_fdm_weights import UpwindSchemeN4FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n6_fdm_weights import UpwindSchemeN6FDMWeights
from methods.fdm.schemes.upwind.upwind_scheme_n8_fdm_weights import UpwindSchemeN8FDMWeights

from methods.fdm.schemes.downwind.downwind_scheme_n2_fdm_weights import DownwindSchemeN2FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n4_fdm_weights import DownwindSchemeN4FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n6_fdm_weights import DownwindSchemeN6FDMWeights
from methods.fdm.schemes.downwind.downwind_scheme_n8_fdm_weights import DownwindSchemeN8FDMWeights


@unique
class FDMEnum(Enum):
    CENTRAL_N2 = (2, CentralSchemeN2FDMWeights())
    CENTRAL_N4 = (4, CentralSchemeN4FDMWeights())
    CENTRAL_N6 = (6, CentralSchemeN6FDMWeights())
    CENTRAL_N8 = (8, CentralSchemeN8FDMWeights())
    UPWIND_N2 = (2, UpwindSchemeN2FDMWeights())
    UPWIND_N4 = (4, UpwindSchemeN4FDMWeights())
    UPWIND_N6 = (6, UpwindSchemeN6FDMWeights())
    UPWIND_N8 = (8, UpwindSchemeN8FDMWeights())
    DOWNWIND_N2 = (2, DownwindSchemeN2FDMWeights())
    DOWNWIND_N4 = (4, DownwindSchemeN4FDMWeights())
    DOWNWIND_N6 = (6, DownwindSchemeN6FDMWeights())
    DOWNWIND_N8 = (8, DownwindSchemeN8FDMWeights())