from typing import List, Tuple, Callable
import matplotlib.pyplot as plt

# номер студента в бригаде
K = 3

# вариант #2: y(x) = x^6 - x
def func(x):
    return x**6 - x

# x \in [-3; 3]
A = -3.
B = 3.


def main():
    # Подготовить исходные данные к заданию
    data = {
        'sample_points': [],
        'tabulated_points': [],
        'L_1': [],
        'L_3': [],
        'L_6': [],
        'nodes_1': [],
        'nodes_3': [],
        'nodes_6': []
    }
    h = (B - A) / 12.
    for i in range(13):
        x_i = A + (h * i)
        data['sample_points'].append((x_i, func(x_i)))

    # Выделить в таблице узлы интерполирования
    for l in [1, 3, 6]:
        data[f'nodes_{l}'] = [data['sample_points'][K+i] for i in range(l+1)]

    # составить таблицу значений соответствующих многочленов
    # порядка не выше первого, третьего, шестого на отрезке с [a, b] с шагом h1
    h1 = (B - A) / 50.
    for i in range(51):
        x_i = A + (h1 * i)
        # Протабулировать заданную функцию на том же отрезке
        data['tabulated_points'].append((x_i, func(x_i)))
        for l in [1, 3, 6]:
            nodes = data[f'nodes_{l}']
            L = get_L(nodes, l)
            data[f'L_{l}'].append((x_i, L(x_i)))

    print_header()
    for i in range(51):
        print_row(data, i)

    for l in [1, 3, 6]:
        show_plot(data, l)


# Интерполяционный многочлен в форме Лагранжа
def get_L(nodes: List[Tuple[float, float]], n: int) -> Callable[[float], float]:
    def internal(x):
        result = 0.
        for i in range(n+1):
            dividend = 1.
            divisor = 1.
            for j in range(n+1):
                if i != j:
                    dividend *= (x - nodes[j][0])
                    divisor *= (nodes[i][0] - nodes[j][0])
            result += nodes[i][1] * (dividend / divisor)
        return result
    return internal


# Построить и напечатать графики функции и интерполяционных многочленов,
# используя полученные значения. На каждом рисунке привести
# график одного интерполяционного многочлена и график заданной функции.
def show_plot(data, l) -> None:
    # plt.ylim([min(p[1] for p in func_points) - 50, max(p[1] for p in func_points) + 50])
    plt.ylim([-5, 15])
    plt.xlim([-5, 5])
    plt.plot(*zip(*data['tabulated_points']), color='red', label='y(x)')
    plt.plot(*zip(*data[f'L_{l}']), color='blue', label=f'L_{l}')
    plt.plot(*zip(*data[f'nodes_{l}']), ls='', marker='s', color='black', markerfacecolor='none', label='y_k')
    plt.legend()
    plt.show()

COL_SIZE = 8
X_DIGITS = 4

def print_header():
    header_template = ' '.join(['{:>{align}}'] * 5)
    print(header_template.format(
        'x', 'y(x)', 'L_1(x)', 'L_3(x)', 'L_6(x)',
        align=COL_SIZE
    ))

def print_row(data, i):
    row_template = ' '.join(['{:>{align}.{digits}f}'] * 5)
    print(row_template.format(
        data['tabulated_points'][i][0],
        data['tabulated_points'][i][1],
        data['L_1'][i][1],
        data['L_3'][i][1],
        data['L_6'][i][1],
        digits=X_DIGITS,
        align=COL_SIZE
    ))

if __name__ == "__main__":
    main()
