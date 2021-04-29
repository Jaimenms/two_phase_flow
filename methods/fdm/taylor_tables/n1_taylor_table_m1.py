from .taylor_table import TaylorTable


class N1TaylorTableM1(TaylorTable):

    @staticmethod
    def lambdas():

        out = [
            [None, None, None],
            [None, None, None],
            [None, None, None],
        ]

        out[0][0] = lambda x: -1 / (-x[0] + x[1])
        out[0][1] = lambda x: 1 / (-x[0] + x[1])
        out[1][0] = lambda x: -1 / (-x[0] + x[1])
        out[1][1] = lambda x: 1 / (-x[0] + x[1])

        return out