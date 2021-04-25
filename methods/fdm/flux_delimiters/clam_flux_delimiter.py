from .flux_delimiter import FluxDelimiter


class ClamFluxDelimiter(FluxDelimiter):

    def __call__(self, phi_p: float, tetha_f: float, tetha_p):

        if 0 < phi_p <= 1:
            return (tetha_f - tetha_p**2)/(tetha_p*(1-tetha_p))*phi_p - (tetha_f - tetha_p)/(tetha_p*(1-tetha_p))*phi_p**2
        else:
            return phi_p
