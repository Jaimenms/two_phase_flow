from ..scheme_fdm_mixin import SchemeFDMMixin


class CentralSchemeFDMMixin(SchemeFDMMixin):

    @staticmethod
    def elements(N: int, L: int, i: int):

        if i < int(N/2):
            l = i
            ini = 0
            fini = N + 1
        elif i >= L-int(N/2):
            l = i - L
            ini = -(N + 1)
            fini = L
        else:
            l = int(N / 2)
            ini = i - int(N / 2)
            fini = i + int(N / 2) + 1

        return l, ini, fini
