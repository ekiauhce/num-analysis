import csv

import matplotlib.pyplot as plt
import numpy as np
from numpy import poly1d


def func(x):
    return x**6 - x


def get_table(a, b):
    step = 10
    h = (b - a) / step
    answer = list()
    if a < b:
        for i in range(step + 1):
            x_i = a + i * h
            answer.append([x_i, func(x_i)])

        return answer
    else:
        print("a>=b")


def calculate(nodes_interpolation, m):
    part_nodes = nodes_interpolation[0:m]
    A = []
    b = []

    for i in range(0, m):
        coefficients = list()
        d = 0
        for j in range(m):
            coefficients.append(sum([pow(x[0], i+d) for x in part_nodes]))
            d += 1

        A.append(coefficients)
        b.append(sum(pow(e[0], i) * e[1] for e in part_nodes))

    x = list(np.linalg.solve(A, b))
    x.reverse()

    return poly1d(x)


def save_plot(poly, interpolation_nodes, n, path, a=-3, b=3):
    h = (b-a) / 50
    x = np.arange(a, b + h, h)
    nodes_x = [node[0] for node in interpolation_nodes]
    nodes_y = [node[1] for node in interpolation_nodes]

    plt.ylim([-50, 800])
    plt.plot(x, func(x), color='red', label='y(x)')
    plt.plot(nodes_x, nodes_y, ls='', color='black', marker='s', markerfacecolor='none', label='y_k')
    plt.plot(x, poly(x), color='blue', label=f'P_{n}')
    plt.legend()
    plt.savefig(f'{path}poly_{n}.png')
    plt.close()


if __name__ == "__main__":
    path = "result/"
    a = -3
    b = 3
    nodes = get_table(a, b)
    stages = [1, 2, 5, 6]

    # записываем узлы интерполирования
    with open(f"{path}nodes.csv", "w+", newline='') as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerow(['x'] + [round(x[0], 3) for x in nodes])
        csvWriter.writerow(['y'] + [round(y[1], 3) for y in nodes])

    data = {}
    for k in stages:
        t = [nodes[i] for i in range(0, k+1)]
        p = calculate(nodes, k+1)
        data[f'poly_{k}'] = p
        data[f'table_{k}'] = t

    h = (b - a) / 50
    comparison = []
    for i in range(51):
        x = a + i * h
        poly_val = []
        for j in stages:
            poly_val.append(data[f'poly_{j}'](x))
        cur = [x, func(x), *poly_val]
        comparison.append([round(el, 3) for el in cur])

    with open(f"{path}comparison.csv", "w+", newline='') as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerow(['x', 'f', 'P_1(x)', 'P_2(x)', 'P_5(x)', 'P_6(x)'])
        csvWriter.writerows(comparison)

    # строим графики для анализа интерполяции
    for i in stages:
        poly = data[f'poly_{i}']
        save_plot(poly, data[f'table_{i}'], i, path)
        print(f"P_{len(poly.coef)-1}(x) =", end=" ")
        for x_power, coef in zip(range(len(poly.coef)-1, -1, -1), poly.coef):
            print(f"{coef}x^{x_power}", end=" ")
        print()

