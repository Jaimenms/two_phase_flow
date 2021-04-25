from .flux_delimiter import FluxDelimiter


class SmartFluxDelimiter(FluxDelimiter):

    def __call__(self, phi_p: float, tetha_f: float, tetha_p):

        if 0 < phi_p <= tetha_p/3:
            return tetha_f/tetha_p*phi_p*(1-3*tetha_p+2*tetha_f)/(1-tetha_p)
        elif tetha_p/3 < phi_p <= tetha_p/tetha_f*(1+tetha_f-tetha_p):
            return (tetha_f/tetha_p)*((1-tetha_f)/(1-tetha_p))*phi_p + (tetha_f/(1-tetha_p))*(tetha_f-tetha_p)
        elif tetha_p/tetha_f*(1+tetha_f-tetha_p) < phi_p <= 1:
            return 1.
        else:
            return phi_p
