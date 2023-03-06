import math
import copy

A = [
    [3.51, 0.17, 3.68, -0.28],
    [4.52, 2.11, 6.63, -0.12],
    [-2.11, 3.17, 1.06, -0.15],
    [3.17, 1.81, -3.17, 0.22]
]

B = [
    [0.75],
    [1.11],
    [0.21],
    [0.05]
]
N = 4 # Размерность квадратной матрицы A

def mult(A, B): # матричное произведение
    n = len(A)
    m = len(A[0]) # or len(B)
    p = len(B[0])

    C = [[0.0 for _ in range(p)] for _ in range(n)] # n * p
    for i in range(n):
        for j in range(p):
            s = 0.0
            for k in range(m):
                s += A[i][k] * B[k][j]
            C[i][j] = s
    return C

def transpose(A): # транспонирование матрицы
    n = len(A)
    m = len(A[0])
    A_t = [[0.0 for _ in range(n)] for _ in range(m)] # m * n
    for i in range(len(A)):
        for j in range(len(A[0])):
            A_t[j][i] = A[i][j]
    return A_t

def delta(A): # расчет коэффициента delta
    C = mult(transpose(A), A)
    n = len(C)
    m = len(C[0])
    return min(
        max(sum(
            abs(C[i][j]) for j in range(m))
            for i in range(n)
        ),
        max(sum(
            abs(C[i][j]) for i in range(m))
            for j in range(n)
        ),
        math.sqrt(sum(
            sum(
                C[i][j] * C[i][j] for j in range(m)
            ) for i in range(n)
        ))
    )

def div_scalar(A, scalar): # деление матрицы на скаляр
    n = len(A)
    m = len(A[0])
    R = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            R[i][j] = A[i][j] / scalar
    return R

def mult_scalar(A, scalar): # умножение матрицы на скаляр
    n = len(A)
    m = len(A[0])
    R = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            R[i][j] = A[i][j] * scalar
    return R

def sub(A, B): # разность двух матриц
    n = len(A)
    m = len(A[0])
    R = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            R[i][j] = A[i][j] - B[i][j]
    return R

def add(A, B): # сумма двух матриц
    n = len(A)
    m = len(A[0])
    R = [[0.0 for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            R[i][j] = A[i][j] + B[i][j]
    return R

def E(n): # единичная матрица размера n
    R = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                R[i][j] = 1
    return R

def residual(X1, X2): # невязка
    return max(abs(X1[i][0] - X2[i][0]) for i in range(len(X1)))

COL_SIZE = 12
X_DIGITS = 8

def print_header():
    header_template = ' '.join(['{:>{align}}'] * N) + \
        ' {:>{align}} {:>{align}}'
    print(header_template.format(
        *[f'X{i+1}' for i in range(N)],
        'alpha',
        'iterations',
        align=COL_SIZE
    ))

def print_row(X, alpha, iterations):
    row_template = ' '.join(['{:>{align}.{digits}f}'] * N) + \
        ' {:>{align}.2f} {:>{align}}'
    print(row_template.format(
        *[x[0] for x in X],
        alpha,
        iterations,
        digits=X_DIGITS,
        align=COL_SIZE
    ))


def main():
    EPSILON = 1e-3

    print_header()
    # alpha от 0.1 до 1.9 с шагом в 0.1
    for alpha in (0.1 * x for x in range(1, 20)):
        G = sub(
            E(N),
            mult_scalar(
                div_scalar(mult(transpose(A), A), delta(A)),
                alpha
            )
        )
        F = mult_scalar(
            div_scalar(mult(transpose(A), B), delta(A)),
            alpha
        )

        X_prev = copy.deepcopy(F)
        X = add(mult(G, X_prev), F)
        iterations = 1
        while residual(X, X_prev) > EPSILON:
            X_prev = X
            X = add(mult(G, X_prev), F)
            iterations += 1
        print_row(X, alpha, iterations)

if __name__ == "__main__":
    main()
