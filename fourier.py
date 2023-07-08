from math import pi, cos, sin, exp
import matplotlib.pyplot as plt


T = 2 * pi
def func(x):
    if -pi <= x <= 0:
        return (x / pi) + 1
    elif 0 < x: # <= pi
        return (-x / pi) + 1

# пример
# T = 4
# def func(x):
#     if -2 <= x <= -1:
#         return 1 / x
#     elif -1 < x <= 0:
#         return x**2
#     elif 0 < x <= 2:
#         return exp(-x)

OMEGA = 2*pi / T

def get_h(n):
    return T / (2*n)

def sample(n):
    points = {}
    h = get_h(n)
    for k in range(-n+1, n+1):
        x_k = (-T/2) + (h * (k+n))
        points[k] = (x_k, func(x_k))

    return points

def a(l, n, points):
    sum = 0.
    h = get_h(n)
    for k in range(-n+1, n+1):
        y_k = points[k][1]
        sum += y_k * cos(l * OMEGA * k * h)
    return (h/T) * sum

def b(l, n, points):
    sum = 0.
    h = get_h(n)
    for k in range(-n+1, n+1):
        y_k = points[k][1]
        sum += y_k * sin(l * OMEGA * k * h)
    return (h/T) * sum

def alpha(l, n, points):
    if l == 0 or l == n:
        return a(l, n, points)
    return 2 * a(l, n, points)

def beta(l, n, points):
    if l == 0 or l == n:
        return 0
    return 2 * b(l, n, points)

def phi(n, points, x):
    sum = 0.
    for l in range(1, n+1):
        sum += \
            alpha(l, n, points) * cos(l * OMEGA * x) + \
            beta(l, n, points) * sin(l * OMEGA * x)
    return alpha(0, n, points) + sum

X_DIGITS = 3

def print_header(cols, min_col_size):
    header_template = ' | '.join(
        [f'{{:>{max(min_col_size, len(c))}}}' for c in cols]
    )
    print(header_template.format(
        *cols
    ))

def print_row(cols, row, min_col_size, digits=None):
    row_template = ' | '.join(
        [f'{{:>{max(min_col_size, len(c))}.{{digits}}f}}' for c in cols]
    )
    print(row_template.format(
        *row,
        digits=digits if digits else X_DIGITS,
    ))

nodes_cols = [
    'x_k',
    'y_k'
]
for n in [4, 8, 16]:
    points = sample(n)
    print(f'n = {n}')
    print_header(nodes_cols, min_col_size=8)
    for p in points.values():
        print_row(
            cols=nodes_cols,
            row=p,
            min_col_size=8,
        )
    print()


coefs_cols = [
    'alpha_l',
    'beta_l'
]
for n in [4, 8, 16]:
    points = sample(n)
    print(f'n = {n}, alpha_0 = {alpha(0, n, points):.3f}')
    col_size = 5
    print_header(coefs_cols, min_col_size=col_size)
    for l in range(1, n+1):
        print_row(
            cols=coefs_cols,
            row=[alpha(l, n, points), beta(l, n, points)],
            min_col_size=col_size,
            digits=3,
        )
    print()

tab_cols = [
    'x',
    'y(x)',
    'phi(x), n=4',
    'phi(x), n=8',
    'phi(x), n=16'
]
TAB_MIN_COL_SIZE = 8

print_header(tab_cols, TAB_MIN_COL_SIZE)

h = T / 50
for i in range(0, 51):
    x_i = (-T/2) + (h * i)
    print_row(
        cols=tab_cols,
        min_col_size=TAB_MIN_COL_SIZE,
        row=[
            x_i,
            func(x_i),
            phi(4, sample(4), x_i),
            phi(8, sample(8), x_i),
            phi(16, sample(16), x_i)
        ]
    )


for n in [4, 8, 16]:
    plt.ylim([-0.25, 1.2])
    plt.xlim([-0.25 + -T/2, T/2 + 0.25])
    x_values = [(-T/2) + ((T/50) *i) for i in range(0, 51)]
    plt.plot(x_values, [func(x) for x in x_values], color='red', label='y(x)')
    plt.plot(
        x_values,
        [phi(n, sample(n), x) for x in x_values],
        color='blue',
        label=f'phi(x), n={n}'
    )
    plt.plot(
        *zip(*sample(n).values()),
        ls='',
        marker='s',
        color='black',
        markerfacecolor='none',
        label='y_k'
    )
    plt.legend()
    plt.show()
