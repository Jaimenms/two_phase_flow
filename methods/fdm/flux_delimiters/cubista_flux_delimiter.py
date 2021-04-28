from .flux_delimiter import FluxDelimiter
import numpy as np


class CubistaFluxDelimiter(FluxDelimiter):

    def __call__(self, phi_p: np.ndarray, tetha_f: np.ndarray, tetha_p: np.ndarray) -> np.ndarray:

        return np.select(
            condlist=[
                phi_p<=3/4*tetha_p,
                phi_p<=(1 + 2*(tetha_f - tetha_p))/(2*tetha_f - tetha_p)*tetha_p,
                phi_p<=1,
                phi_p>1,
            ],
            choicelist=[
                (1+(tetha_f-tetha_p)/(2*(1-tetha_p)))*tetha_f/tetha_p*phi_p,
                tetha_f/tetha_p*(1-tetha_f)/(1-tetha_p)*phi_p+tetha_f*(tetha_f-tetha_p)/(1-tetha_p),
                1-(1-tetha_f)/(2*(1-tetha_p))*(1-phi_p),
                phi_p
            ],
            default=0.
        )