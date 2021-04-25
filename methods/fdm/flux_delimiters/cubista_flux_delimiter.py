from .flux_delimiter import FluxDelimiter


class CubistaFluxDelimiter(FluxDelimiter):

    def __call__(self, phi_p: float, tetha_f: float, tetha_p):

        if 0 < phi_p <= 3/4*tetha_p:

            return (1+(tetha_f-tetha_p)/(2*(1-tetha_p)))*tetha_f/tetha_p*phi_p

        elif 3/4*tetha_p < phi_p <= (1 + 2*(tetha_f - tetha_p))/(2*tetha_f - tetha_p)*tetha_p:

            return tetha_f/tetha_p*(1-tetha_f)/(1-tetha_p)*phi_p+tetha_f*(tetha_f-tetha_p)/(1-tetha_p)

        elif (1 + 2*(tetha_f - tetha_p))/(2*tetha_f - tetha_p)*tetha_p < phi_p <= 1:

            return 1-(1-tetha_f)/(2*(1-tetha_p))*(1-phi_p)

        else:

            return phi_p
