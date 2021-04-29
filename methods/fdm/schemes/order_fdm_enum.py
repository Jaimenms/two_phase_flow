from enum import Enum, unique

from methods.fdm.schemes.scheme_m1_fdm_enum import SchemeM1FDMEnum
from methods.fdm.schemes.scheme_m2_fdm_enum import SchemeM2FDMEnum


@unique
class OrderFDMEnum(Enum):
    FIRST = (1, SchemeM1FDMEnum)
    SECOND = (2, SchemeM2FDMEnum)
