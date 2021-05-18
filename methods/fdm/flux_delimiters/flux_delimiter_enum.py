from enum import Enum, unique
from .minmod_flux_delimiter import MinModFluxDelimiter
from .minmod2_flux_delimiter import MinMod2FluxDelimiter
from .smart_flux_delimiter import SmartFluxDelimiter
from .smart2_flux_delimiter import Smart2FluxDelimiter
from .clam_flux_delimiter import ClamFluxDelimiter
from .cubista_flux_delimiter import CubistaFluxDelimiter
from .cubista2_flux_delimiter import Cubista2FluxDelimiter


@unique
class FluxDelimiterEnum(Enum):
    MINMOD = MinModFluxDelimiter()
    MINMOD2 = MinMod2FluxDelimiter()
    SMART = SmartFluxDelimiter()
    SMART2 = Smart2FluxDelimiter()
    CLAM = ClamFluxDelimiter()
    CUBISTA = CubistaFluxDelimiter()
    CUBISTA2 = Cubista2FluxDelimiter()
