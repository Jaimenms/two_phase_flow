from .flux_delimiter import FluxDelimiter
import numpy as np
from numba import vectorize, float32, float64


class MinMod2FluxDelimiter(FluxDelimiter):

    @staticmethod
    @vectorize([float32(float32, float32,float32),
                float64(float64, float64,float64)])
    def __call__(phi_p, tetha_f, tetha_p):
        if phi_p<=tetha_p:
            return tetha_f/tetha_p*phi_p
        elif phi_p<=1:
            return ((1-tetha_f)*phi_p + (tetha_f-tetha_p))/(1-tetha_p)
        else:
            return phi_p
