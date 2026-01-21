import sys

import numpy as np
import time
from sympy import lambdify, symbols

from app.Functions_cos import wn_cos, right_part_eq1_cos, eq1_integral_cos

a = 5  # m
b = 6  # m
E = 27500  # MPa
mu = 0.2
h = 0.2  # m
q = 0.01  # MPa

D = E * h ** 3 / (12 * (1 - mu ** 2))

n_approx = 6  # Номер приближения
n_wi = 1  # Кол-во неизвестных при фиксированных i и j
indexes = 2  # Всего 2 индекса: i и j

np.set_printoptions(threshold=sys.maxsize)  # чтобы выводилась полностью вся матрица

t1 = time.time()

# Сбор основной системы уравнений
number_variables = n_wi * (n_approx + 1) ** indexes  # wi
left_part: np.ndarray = np.zeros((number_variables, number_variables))
right_part: np.ndarray = np.zeros(number_variables)

equation_position = 0

for k in range(n_approx + 1):
    for l in range(n_approx + 1):
        # Сначала индексы k, l по индексам в уравнении
        position_wij_in_eq1 = 0

        right_part[0 + n_wi * equation_position] = right_part_eq1_cos(k, l, a, b, q)[0]

        print('k = ', k)
        print('l = ', l)

        for i in range(n_approx + 1):
            for j in range(n_approx + 1):
                # print('i = ', i)
                # print('j = ', j)
                w1_in_Wn = wn_cos(i, j, a, b)

                left_part[0 + n_wi * equation_position][position_wij_in_eq1] = eq1_integral_cos(k, l, a, b, D, w1_in_Wn)[0]
                position_wij_in_eq1 += 1
                # print(left_part)
                # print(right_part)

        equation_position += 1

print(left_part)
print(right_part)

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
        w1_in_Wn = wn_cos(i, j, a, b)
        # print('n_wi * equation_position = ', n_wi * equation_position)
        Wn_raw = ((w1_in_Wn) * vector_w[0 + n_wi * equation_position])

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

