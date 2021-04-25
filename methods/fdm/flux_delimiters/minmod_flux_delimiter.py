from .flux_delimiter import FluxDelimiter


class MinModFluxDelimiter(FluxDelimiter):

    def __call__(self, phi_p: float, tetha_f: float, tetha_p):

        if 0 < phi_p <= tetha_p:
            return tetha_f/tetha_p*phi_p
        elif tetha_p < phi_p <= 1:
            return ((1-tetha_f)*phi_p + (tetha_f-tetha_p))/(1-tetha_p)
        else:
            return phi_p
