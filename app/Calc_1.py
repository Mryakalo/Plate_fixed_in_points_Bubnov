import sys
from queue import Queue
from typing import Callable

import numpy as np
import time

from numpy.f2py.cfuncs import includes
from sympy import lambdify, symbols, Expr
from multiprocessing import Manager, Pool

from app.Functions import wn, eq1_integral, right_part_eq1, eq2_integral, right_part_eq2, eq3_integral, right_part_eq3, \
    eq4_integral, right_part_eq4
from app.Plot_diagram import plot_graph

np.set_printoptions(threshold=sys.maxsize)  # чтобы выводилась полностью вся матрица

a = 10  # m
b = 10  # m
E = 27500  # MPa
mu = 0.2
h = 0.2  # m
q = 0.01  # MPa

D = E * h ** 3 / (12 * (1 - mu ** 2))

n_approx = 3  # Номер приближения
n_wi = 4  # Кол-во неизвестных при фиксированных i и j
indexes = 2  # Всего 2 индекса: i и j

NUMBER_OF_QUEUES = (n_approx + 2) ** indexes  # число очередей для сохранения результата из процессов


def store_integral(
    row,
    column: int,
    k: int,
    l: int,
    a: int,
    b: int,
    D: float,
    w_in_Wn: Expr,
    integral_function: Callable,
    queues: list[Queue],
) -> None:

    value = integral_function(k, l, a, b, D, w_in_Wn)

    # выбираем одну очередь из списка очередей по номеру ряда, чтобы очереди не забивались
    queue = queues[row % NUMBER_OF_QUEUES]
    queue.put_nowait((row, column, value))


def main() -> None:
    print(f'\nApproximation: {n_approx}\n')
    equation_position = 0
    t1 = time.time()

    # Сбор основной системы уравнений
    number_variables = n_wi * (n_approx + 1) ** indexes  # wi
    left_part: np.ndarray = np.zeros((number_variables, number_variables))
    right_part: np.ndarray = np.zeros(number_variables)

    with Manager() as manager:
        queues = [manager.Queue() for _ in range(NUMBER_OF_QUEUES)]

        for k in range(n_approx + 1):
            print('k = ', k)
            for l in range(n_approx + 1):
                # Сначала индексы k, l по индексам в уравнении
                position_wij_in_eq1 = 0
                position_wij_in_eq2 = 0
                position_wij_in_eq3 = 0
                position_wij_in_eq4 = 0
                right_part[0 + n_wi * equation_position] = right_part_eq1(k, l, a, b, q)[0]
                right_part[1 + n_wi * equation_position] = right_part_eq2(k, l, a, b, q)[0]
                right_part[2 + n_wi * equation_position] = right_part_eq3(k, l, a, b, q)[0]
                right_part[3 + n_wi * equation_position] = right_part_eq4(k, l, a, b, q)[0]
                print('l = ', l)

                for i in range(n_approx + 1):
                    # print('i = ', i)
                    for j in range(n_approx + 1):
                        # print('j = ', j)
                        w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn = wn(i, j, a, b)

                        with Pool(16) as pool:  # подсчёт данных через 16 процессов и складывание их в разные очереди
                            pool.starmap(
                                store_integral,
                                [
                                    # порядок аргументов важен и должен быть такой же как и в методе store_integral
                                    (0 + n_wi * equation_position, position_wij_in_eq1, k, l, a, b, D, w1_in_Wn, eq1_integral, queues),
                                    (0 + n_wi * equation_position, position_wij_in_eq1 + 1, k, l, a, b, D, w2_in_Wn, eq1_integral, queues),
                                    (0 + n_wi * equation_position, position_wij_in_eq1 + 2, k, l, a, b, D, w3_in_Wn, eq1_integral, queues),
                                    (0 + n_wi * equation_position, position_wij_in_eq1 + 3, k, l, a, b, D, w4_in_Wn, eq1_integral, queues),
                                    (1 + n_wi * equation_position, position_wij_in_eq2, k, l, a, b, D, w1_in_Wn, eq2_integral, queues),
                                    (1 + n_wi * equation_position, position_wij_in_eq2 + 1, k, l, a, b, D, w2_in_Wn, eq2_integral, queues),
                                    (1 + n_wi * equation_position, position_wij_in_eq2 + 2, k, l, a, b, D, w3_in_Wn, eq2_integral, queues),
                                    (1 + n_wi * equation_position, position_wij_in_eq2 + 3, k, l, a, b, D, w4_in_Wn, eq2_integral, queues),
                                    (2 + n_wi * equation_position, position_wij_in_eq3, k, l, a, b, D, w1_in_Wn, eq3_integral, queues),
                                    (2 + n_wi * equation_position, position_wij_in_eq3 + 1, k, l, a, b, D, w2_in_Wn, eq3_integral, queues),
                                    (2 + n_wi * equation_position, position_wij_in_eq3 + 2, k, l, a, b, D, w3_in_Wn, eq3_integral, queues),
                                    (2 + n_wi * equation_position, position_wij_in_eq3 + 3, k, l, a, b, D, w4_in_Wn, eq3_integral, queues),
                                    (3 + n_wi * equation_position, position_wij_in_eq4, k, l, a, b, D, w1_in_Wn, eq4_integral, queues),
                                    (3 + n_wi * equation_position, position_wij_in_eq4 + 1, k, l, a, b, D, w2_in_Wn, eq4_integral, queues),
                                    (3 + n_wi * equation_position, position_wij_in_eq4 + 2, k, l, a, b, D, w3_in_Wn, eq4_integral, queues),
                                    (3 + n_wi * equation_position, position_wij_in_eq4 + 3, k, l, a, b, D, w4_in_Wn, eq4_integral, queues),
                                ]
                            )
                        position_wij_in_eq1 += 4
                        position_wij_in_eq2 += 4
                        position_wij_in_eq3 += 4
                        position_wij_in_eq4 += 4

                equation_position += 1

            t_temp = time.time()
            time_spend = int(t_temp - t1)
            total = int(time_spend / (k + 1) * (n_approx + 1))
            eta = total - time_spend

            print(f'time spend: {time_spend} seconds')
            print(f'time to finish: {eta} seconds')
            print(f'expected total time: {total} seconds')

        # Сборка left_part. Данные в очередях лежат в случайном порядке, в зависимости какой процесс отработал быстрее.
        for queue in queues:
            while not queue.empty():
                element = queue.get_nowait()
                row = element[0]
                column = element[1]
                result = element[2]
                left_part[row][column] = result

    # print(left_part)
    # print(right_part)
    vector_w: np.ndarray = np.linalg.solve(left_part, right_part)
    print(vector_w)

    x, y = symbols('x y')
    Wn = 0
    equation_position = 0
    sign = 1

    for i in range(n_approx + 1):
        for j in range(n_approx + 1):
            print('i = ', i)
            print('j = ', j)
            w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn = wn(i, j, a, b)
            # print('n_wi * equation_position = ', n_wi * equation_position)
            Wn_raw = ((w1_in_Wn) * vector_w[0 + n_wi * equation_position] +
                      (w2_in_Wn) * vector_w[1 + n_wi * equation_position] +
                      (w3_in_Wn) * vector_w[2 + n_wi * equation_position] +
                      (w4_in_Wn) * vector_w[3 + n_wi * equation_position])

            # print('w1_in_Wn = ', w1_in_Wn)
            # print('w2_in_Wn = ', w2_in_Wn)
            # print('w3_in_Wn = ', w3_in_Wn)
            # print('w4_in_Wn = ', w4_in_Wn)

            Wn = (Wn + Wn_raw)

            Wn_lambda = lambdify([x, y], Wn_raw)
            W_middle = Wn_lambda(
                a / 2,
                b / 2
            )
            print('Часть от перемещения, мм', W_middle * 1000)
            equation_position += 1

    Wn_lambda = lambdify([x, y], Wn)
    W_middle = Wn_lambda(
        a / 2,
        b / 2
    )
    print('Перемещения в центре, мм', W_middle * 1000)
    W_edge = Wn_lambda(
        0,
        b / 2
    )
    print('Перемещения на краю, мм', W_edge * 1000)
    W_edge = Wn_lambda(
        a / 2,
        0
    )
    print('Перемещения на краю, мм', W_edge * 1000)
    W_edge = Wn_lambda(
        a,
        b
    )
    print('Перемещения на краю, мм', W_edge * 1000)

    t2 = time.time()
    t = t2 - t1
    print('time', t)

    plot_graph(Wn_lambda, a, b)


if __name__ == "__main__":
    main()