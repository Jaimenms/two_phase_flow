from ..scheme_fdm_mixin import SchemeFDMMixin

class UpwindSchemeFDMMixin(SchemeFDMMixin):

    @staticmethod
    def elements(N: int, L: int, i: int):

        if i < L - N:
            l = 0
            ini = i
            fini = i + N +1
        else:
            l = -(L - i)
            ini = -(N + 1)
            fini = L

        return l, ini, fini
