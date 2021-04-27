from sympy import Symbol


def fornberg(M, N, alpha, x0):
    """
    Fornberg Algorithm
    Fornberg B. (1988) Generation of Finite Difference Formulas on Arbitrarily Spaced Grids,
    Mathematics of Computation 51, no. 184 : 699-706.
    https://www.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0935077-0/S0025-5718-1988-0935077-0.pdf
    :param M: degree of maximum differentiation
    :param N: approximation order
    :param alpha: list of symbols for the points
    :param x0: symbol of the approximation point
    :return:
    """
    delta = []
    for i in range(M + 1):
        delta_n = []
        for n in range(N + 1):
            delta_k = []
            for k in range(N + 1):
                delta_k.append(0)
            delta_n.append(delta_k)
        delta.append(delta_n)
        delta[0][0][0] = 1

    c1 = 1

    for n in range(1, N + 1):
        c2 = 1
        for v in range(0, n):
            c3 = alpha[n] - alpha[v]
            c2 *= c3
            if n <= M:
                delta[n][n - 1][v] = 0
            for m in range(0, min(n, M) + 1):
                delta[m][n][v] = ((alpha[n] - x0) * delta[m][n - 1][v] - m * delta[m - 1][n - 1][v]) / c3

        for m in range(0, min(n, M) + 1):
            delta[m][n][n] = c1 / c2 * (
                        m * delta[m - 1][n - 1][n - 1] - (alpha[n - 1] - x0) * delta[m][n - 1][n - 1])

        c1 = c2

    return delta[M][-1]


def print_taylor_table(M, N, kind="homogeneous"):
    """
    Prints the Taylor Matrix represented as lambda functions to be copy/pasted to the corresponting py file
    :param M: degree of differentiation
    :param N: approximation order
    :param kind: "homogeneous" if all points are equally spaced, "heterogeneous" if not
    :return: NONE
    """
    alpha = []
    for ialpha in range(0, N + 1):
        if kind == "homogeneous":
            alpha.append(Symbol('h') * ialpha)
        else:
            alpha.append(Symbol('x[{}]'.format(ialpha)))

    j = 0
    for jalpha in range(0, N + 1):

        if kind == "homogeneous":
            x0 = Symbol('h') * jalpha
        else:
            x0 = Symbol("x[{}]".format(jalpha))

        out = fornberg(M, N, alpha, x0)

        print("#")
        print("#df/dx({})".format(x0))
        for node, outi in enumerate(out):
            print("out[{}][{}] = lambda x : {}".format(j, node, outi))

        j += 1