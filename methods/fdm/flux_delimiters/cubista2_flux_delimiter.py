from .flux_delimiter import FluxDelimiter
from numba import vectorize, float32, float64


class Cubista2FluxDelimiter(FluxDelimiter):

    @staticmethod
    @vectorize([float32(float32, float32,float32),
                float64(float64, float64,float64)])
    def __call__(phi_p, tetha_f, tetha_p):
        if phi_p<=3/4*tetha_p:
            return (1+(tetha_f-tetha_p)/(2*(1-tetha_p)))*tetha_f/tetha_p*phi_p
        elif phi_p<=(1 + 2*(tetha_f - tetha_p))/(2*tetha_f - tetha_p)*tetha_p:
            return tetha_f/tetha_p*(1-tetha_f)/(1-tetha_p)*phi_p+tetha_f*(tetha_f-tetha_p)/(1-tetha_p)
        elif phi_p<=1:
            return 1-(1-tetha_f)/(2*(1-tetha_p))*(1-phi_p)
        else:
            return phi_p
