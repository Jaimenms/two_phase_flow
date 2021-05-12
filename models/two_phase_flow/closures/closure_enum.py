from enum import Enum, unique
from models.two_phase_flow.closures.annular_closure import AnnularClosure
from models.two_phase_flow.closures.bubble_closure import BubbleClosure
from models.two_phase_flow.closures.intermittent_closure import IntermittentClosure
from models.two_phase_flow.closures.stratified_closure import StratifiedClosure
from models.two_phase_flow.closures.automatic_closure import AutomaticClosure


@unique
class ClosureEnum(Enum):
    ANNULAR = AnnularClosure.model_closure
    BUBBLE = BubbleClosure.model_closure
    INTERMITTENT = IntermittentClosure.model_closure
    STRATIFIED = StratifiedClosure.model_closure
    AUTOMATIC = AutomaticClosure.model_closure