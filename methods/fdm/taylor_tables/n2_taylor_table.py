from .taylor_table import TaylorTable


class N2TaylorTable(TaylorTable):

    @staticmethod
    def lambdas():

        out = [
            [None, None, None],
            [None, None, None],
            [None, None, None],
        ]

        out[0][0] = lambda x: (-1 - (-x[0] + x[2])/(-x[0] + x[1]))/(-x[0] + x[2])
        out[0][1] = lambda x: (-x[0] + x[2])/((-x[0] + x[1])*(-x[1] + x[2]))
        out[0][2] = lambda x: -(-x[0] + x[1])/((-x[0] + x[2])*(-x[1] + x[2]))
        out[1][0] = lambda x: -(-x[1] + x[2])/((-x[0] + x[1])*(-x[0] + x[2]))
        out[1][1] = lambda x: (-1 + (-x[1] + x[2])/(-x[0] + x[1]))/(-x[1] + x[2])
        out[1][2] = lambda x: (-x[0] + x[1])/((-x [0] + x[2])*(-x[1] + x[2]))
        out[2][0] = lambda x: -(x[1] - x[2]) / ((-x[0] + x[1]) * (-x[0] + x[2]))
        out[2][1] = lambda x: -(-x[0] + x[2]) / ((-x[0] + x[1]) * (-x[1] + x[2]))
        out[2][2] = lambda x: (-x[0] + x[1]) * ((-x[0] + x[2]) / (-x[0] + x[1]) - (x[1] - x[2]) / (-x[0] + x[1])) / ((-x[0] + x[2]) * (-x[1] + x[2]))

        return out