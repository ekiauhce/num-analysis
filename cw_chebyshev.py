import numpy as np
import math

import numpy as np

A = np.array([
    [ 2, -1,  0],
    [-1,  2, -1],
    [ 0, -1,  2]
], dtype=float)
B = np.array([1, 4, -3], dtype=float)
x0 = np.array([0, 0, 0], dtype=float)

epsilon = 1e-6

lambda_min = 2 - math.sqrt(2)
lambda_max = 2 + math.sqrt(2)

tau_0 = 2 / (lambda_min + lambda_max)
ksi = lambda_min / lambda_max
rho_0 = (1 - ksi) / (1 + ksi)


def n_0(eps):
    return math.log(2 / eps) / (2 * math.sqrt(ksi))


n = round(n_0(epsilon))
print(f'Матрица A = \n{A}')
print(f'Вектор B = \n{B}')

print()


def t(k):
    return (math.cos(2*k - 1) * math.pi) / (2 * n)

def tau(k):
    return tau_0 / (1 + rho_0 * t(k))

def norm(vec):
    s = 0.
    for x_i in vec:
        s += abs(x_i)**2
    return s


rho_1 =  (1 - math.sqrt(ksi)) / (1 + math.sqrt(ksi))

def q(n):
    return (2 * rho_1**n) / (1 + rho_1**(2*n))

x = x0
x_next = None
k = 1
diff_norm = None
print(f'Заданная погрешность epsilon = {epsilon}')
print(f'Оценка количества итераций n = {n}')
while diff_norm is None or diff_norm > epsilon:
    x_next = (B - A.dot(x)) * tau(k) + x
    diff_norm = norm(x_next - x)
    print(f'i = {k}, x = {x}, норма = {diff_norm:.6f}')
    x = x_next
    k += 1

print(f'Ответ: вектор x = {x}')
print(f'невязка = {A.dot(x) - B}')
print()


n = 30
epsilon = q(n)
x = x0
x_next = None
k = 1
diff_norm = None
print(f'Заданное количество итераций n = {n}')
print(f'Оценка погрешности epsilon = {epsilon}')
for k in range(1, n+1):
    x_next = (B - A.dot(x)) * tau(k) + x
    print(f'i = {k}, x = {x}, норма = {norm(x_next - x):.10f}')
    x = x_next

print()
print(f'Ответ: вектор x = {x}')
print()

def seidel(A, b, eps):
    n = len(A)
    x = np.zeros(n)

    converge = False
    k = 0
    while not converge:
        k += 1
        x_new = np.copy(x)
        for i in range(n):
            s1 = sum(A[i][j] * x_new[j] for j in range(i))
            s2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_new[i] = (b[i] - s1 - s2) / A[i][i]

        converge = norm(x_new - x) <= eps
        print(f'i = {k}, x = {x}, норма = {norm(x_new - x):.8f}')
        x = x_new
    print(f'невязка = {A.dot(x) - b}')


seidel(A, B, 1e-6)