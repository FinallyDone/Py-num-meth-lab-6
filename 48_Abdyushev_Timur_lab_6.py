import numpy as np
import math

Var = 1
M = 10
N = 10
ALPHA = 0.5 + 0.1 * Var


def phi_0(t):
    return 1 / (ALPHA * t + 1)


def phi_1(t):
    return 1 / (ALPHA * t + 2)


def alpha(x):
    return 1 / (x + 1)


def d2_alpha(x):
    return -(-2 * x - 2) / (x + 1)**4


def beta(x):
    return -ALPHA / (1 + x)**2


def func(x, t):
    return 2 * (ALPHA**2 - 1) / ((x + ALPHA * t + 1)**3)


def get_Yi0(x):
    return alpha(x)


def get_Yi1(x):
    return alpha(x) + Mu * beta(x) + Mu**2 / 2 * (a**2 * d2_alpha(x) + func(x, 0))


def get_Yij(i, j):
    return Sij[i][j] * matrix_Y[i + 1][j] + 2 * (1 - Sij[i][j]) * matrix_Y[i][j] + \
           Sij[i][j] * matrix_Y[i - 1][j] - matrix_Y[i][j - 1] + Mu**2 * func(X[i], t[j])


if __name__ == '__main__':
    a = 1
    l = 1
    T = 1
    h = l / N
    Mu = T / N
    X = [i * h for i in range(N+1)]
    t = [i * Mu for i in range(N+1)]
    Sij = np.ones((N+1, N+1))
    Sij *= Mu**2 / h**2 * a
    matrix_Y = np.zeros((N+1, M+1))

    for j in range(M+1):
        matrix_Y[0][j] = phi_0(t[j])
        matrix_Y[M][j] = phi_1(t[j])

    for i in range(1, N):
        matrix_Y[i][0] = get_Yi0(X[i])
        matrix_Y[i][1] = get_Yi1(X[i])

    for j in range(1, M):
        for i in range(1, N):
            matrix_Y[i][j+1] = get_Yij(i, j)

    print("\nМатрица решений Y(x):\n")
    for i in range(M, -1, -1):
        print('t = ', i/10, end = "  |")
        for j in range(M+1):
            print("{0:8.4f}".format(matrix_Y[j][i]), end = " ")
        print('')

