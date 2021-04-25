from ..scheme_fdm_mixin import SchemeFDMMixin


class DownwindSchemeFDMMixin(SchemeFDMMixin):

    @staticmethod
    def elements(N: int, L: int, i: int):

        if i > N:
            l = N
            ini = i - N
            fini = i + 1
        else:
            l = i
            ini = 0
            fini = N + 1


        return l, ini, fini
