from enum import Enum, unique
from .minmod_flux_delimiter import MinModFluxDelimiter
from .smart_flux_delimiter import SmartFluxDelimiter
from .clam_flux_delimiter import ClamFluxDelimiter
from .cubista_flux_delimiter import CubistaFluxDelimiter


@unique
class FluxDelimiterEnum(Enum):
    MINMOD = MinModFluxDelimiter()
    SMART = SmartFluxDelimiter()
    CLAM = ClamFluxDelimiter()
    CUBISTA = CubistaFluxDelimiter()
