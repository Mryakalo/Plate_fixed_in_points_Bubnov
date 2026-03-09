import sys

import numpy as np
import time
from sympy import lambdify, symbols, Expr

from app.Functions import wn, eq1_integral, right_part_eq1, eq2_integral, right_part_eq2, eq3_integral, right_part_eq3, \
    eq4_integral, right_part_eq4
from app.Plot_diagram import plot_graph


a = 6  # m
b = 6  # m
E = 27500  # MPa
mu = 0.2
h = 0.2  # m
q = 0.01  # MPa

D = E * h ** 3 / (12 * (1 - mu ** 2))

n_approx = 4  # Номер приближения
n_wi = 4  # Кол-во неизвестных при фиксированных i и j
indexes = 2  # Всего 2 индекса: i и j

np.set_printoptions(threshold=sys.maxsize)  # чтобы выводилась полностью вся матрица

t1 = time.time()
# Сбор основной системы уравнений
number_variables = n_wi * (n_approx + 1) ** indexes  # wi
left_part: np.ndarray = np.zeros((number_variables, number_variables))
right_part: np.ndarray = np.zeros(number_variables)

with open('files_for_left_part/eq1.txt', 'w') as f1:
    print('clear file eq1.txt')
with open('files_for_left_part/eq2.txt', 'w') as f2:
    print('clear file eq2.txt')
with open('files_for_left_part/eq3.txt', 'w') as f3:
    print('clear file eq3.txt')
with open('files_for_left_part/eq4.txt', 'w') as f4:
    print('clear file eq4.txt')


def calc_enq1_integral(
        k_in: int,
        l_in: int,
        a_in: float,
        b_in: float,
        D_in: float,
        w_in_Wn: list[Expr],
) -> None:

    with open('files_for_left_part/eq1.txt', 'a') as file1_for_left_part:

        coef_for_w1 = eq1_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[0])
        file1_for_left_part.write(str(coef_for_w1) + '\n')
        coef_for_w2 = eq1_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[1])
        file1_for_left_part.write(str(coef_for_w2) + '\n')
        coef_for_w3 = eq1_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[2])
        file1_for_left_part.write(str(coef_for_w3) + '\n')
        coef_for_w4 = eq1_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[3])
        file1_for_left_part.write(str(coef_for_w4) + '\n')



def calc_enq2_integral(
        k_in: int,
        l_in: int,
        a_in: float,
        b_in: float,
        D_in: float,
        w_in_Wn: list[Expr],
) -> None:

    with open('files_for_left_part/eq2.txt', 'a') as file2_for_left_part:

        coef_for_w1 = eq2_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[0])
        file2_for_left_part.write(str(coef_for_w1) + '\n')
        coef_for_w2 = eq2_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[1])
        file2_for_left_part.write(str(coef_for_w2) + '\n')
        coef_for_w3 = eq2_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[2])
        file2_for_left_part.write(str(coef_for_w3) + '\n')
        coef_for_w4 = eq2_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[3])
        file2_for_left_part.write(str(coef_for_w4) + '\n')


def calc_enq3_integral(
        k_in: int,
        l_in: int,
        a_in: float,
        b_in: float,
        D_in: float,
        w_in_Wn: list[Expr],
) -> None:

    with open('files_for_left_part/eq3.txt', 'a') as file3_for_left_part:

        coef_for_w1 = eq3_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[0])
        file3_for_left_part.write(str(coef_for_w1) + '\n')
        coef_for_w2 = eq3_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[1])
        file3_for_left_part.write(str(coef_for_w2) + '\n')
        coef_for_w3 = eq3_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[2])
        file3_for_left_part.write(str(coef_for_w3) + '\n')
        coef_for_w4 = eq3_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[3])
        file3_for_left_part.write(str(coef_for_w4) + '\n')


def calc_enq4_integral(
        k_in: int,
        l_in: int,
        a_in: float,
        b_in: float,
        D_in: float,
        w_in_Wn: list[Expr],
) -> None:

    with open('files_for_left_part/eq4.txt', 'a') as file4_for_left_part:

        coef_for_w1 = eq4_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[0])
        file4_for_left_part.write(str(coef_for_w1) + '\n')
        coef_for_w2 = eq4_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[1])
        file4_for_left_part.write(str(coef_for_w2) + '\n')
        coef_for_w3 = eq4_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[2])
        file4_for_left_part.write(str(coef_for_w3) + '\n')
        coef_for_w4 = eq4_integral(k_in, l_in, a_in, b_in, D_in, w_in_Wn[3])
        file4_for_left_part.write(str(coef_for_w4) + '\n')


equation_position = 0
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

                calc_enq1_integral(
                    k_in=k,
                    l_in=l,
                    a_in=a,
                    b_in=b,
                    D_in=D,
                    w_in_Wn=[w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn]
                )

                calc_enq2_integral(
                    k_in=k,
                    l_in=l,
                    a_in=a,
                    b_in=b,
                    D_in=D,
                    w_in_Wn=[w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn]
                )

                calc_enq3_integral(
                    k_in=k,
                    l_in=l,
                    a_in=a,
                    b_in=b,
                    D_in=D,
                    w_in_Wn=[w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn]
                )

                calc_enq4_integral(
                    k_in=k,
                    l_in=l,
                    a_in=a,
                    b_in=b,
                    D_in=D,
                    w_in_Wn=[w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn]
                )

        equation_position += 1


with open('files_for_left_part/eq1.txt', 'r') as file_for_left_part_eq1:
    lines_eq1 = file_for_left_part_eq1.readlines()
with open('files_for_left_part/eq2.txt', 'r') as file_for_left_part_eq2:
    lines_eq2 = file_for_left_part_eq2.readlines()
with open('files_for_left_part/eq3.txt', 'r') as file_for_left_part_eq3:
    lines_eq3 = file_for_left_part_eq3.readlines()
with open('files_for_left_part/eq4.txt', 'r') as file_for_left_part_eq4:
    lines_eq4 = file_for_left_part_eq4.readlines()

equation_position = 0
row_number_file_eq1 = 0
row_number_file_eq2 = 0
row_number_file_eq3 = 0
row_number_file_eq4 = 0
for k in range(n_approx + 1):
    print('k = ', k)
    for l in range(n_approx + 1):
        # Сначала индексы k, l по индексам в уравнении
        position_wij_in_eq1 = 0
        position_wij_in_eq2 = 0
        position_wij_in_eq3 = 0
        position_wij_in_eq4 = 0
        print('l = ', l)

        for i in range(n_approx + 1):
            # print('i = ', i)
            for j in range(n_approx + 1):
                # print('j = ', j)

                left_part[0 + n_wi * equation_position][position_wij_in_eq1] = lines_eq1[row_number_file_eq1]
                row_number_file_eq1 += 1
                position_wij_in_eq1 += 1
                left_part[0 + n_wi * equation_position][position_wij_in_eq1] = lines_eq1[row_number_file_eq1]
                row_number_file_eq1 += 1
                position_wij_in_eq1 += 1
                left_part[0 + n_wi * equation_position][position_wij_in_eq1] = lines_eq1[row_number_file_eq1]
                row_number_file_eq1 += 1
                position_wij_in_eq1 += 1
                left_part[0 + n_wi * equation_position][position_wij_in_eq1] = lines_eq1[row_number_file_eq1]
                row_number_file_eq1 += 1
                position_wij_in_eq1 += 1
                # print(left_part)
                # print(right_part)

                left_part[1 + n_wi * equation_position][position_wij_in_eq2] = lines_eq2[row_number_file_eq2]
                row_number_file_eq2 += 1
                position_wij_in_eq2 += 1
                left_part[1 + n_wi * equation_position][position_wij_in_eq2] = lines_eq2[row_number_file_eq2]
                row_number_file_eq2 += 1
                position_wij_in_eq2 += 1
                left_part[1 + n_wi * equation_position][position_wij_in_eq2] = lines_eq2[row_number_file_eq2]
                row_number_file_eq2 += 1
                position_wij_in_eq2 += 1
                left_part[1 + n_wi * equation_position][position_wij_in_eq2] = lines_eq2[row_number_file_eq2]
                row_number_file_eq2 += 1
                position_wij_in_eq2 += 1
                # print(left_part)
                # print(right_part)

                left_part[2 + n_wi * equation_position][position_wij_in_eq3] = lines_eq3[row_number_file_eq3]
                row_number_file_eq3 += 1
                position_wij_in_eq3 += 1
                left_part[2 + n_wi * equation_position][position_wij_in_eq3] = lines_eq3[row_number_file_eq3]
                row_number_file_eq3 += 1
                position_wij_in_eq3 += 1
                left_part[2 + n_wi * equation_position][position_wij_in_eq3] = lines_eq3[row_number_file_eq3]
                row_number_file_eq3 += 1
                position_wij_in_eq3 += 1
                left_part[2 + n_wi * equation_position][position_wij_in_eq3] = lines_eq3[row_number_file_eq3]
                row_number_file_eq3 += 1
                position_wij_in_eq3 += 1
                # print(left_part)
                # print(right_part)

                left_part[3 + n_wi * equation_position][position_wij_in_eq4] = lines_eq4[row_number_file_eq4]
                row_number_file_eq4 += 1
                position_wij_in_eq4 += 1
                left_part[3 + n_wi * equation_position][position_wij_in_eq4] = lines_eq4[row_number_file_eq4]
                row_number_file_eq4 += 1
                position_wij_in_eq4 += 1
                left_part[3 + n_wi * equation_position][position_wij_in_eq4] = lines_eq4[row_number_file_eq4]
                row_number_file_eq4 += 1
                position_wij_in_eq4 += 1
                left_part[3 + n_wi * equation_position][position_wij_in_eq4] = lines_eq4[row_number_file_eq4]
                row_number_file_eq4 += 1
                position_wij_in_eq4 += 1
                # print(left_part)
                # print(right_part)

        equation_position += 1

# print(left_part)
# print(right_part)

vector_w: np.ndarray = np.linalg.solve(left_part, right_part)
# print(vector_w)

x, y = symbols('x y')
Wn = 0
equation_position = 0
sign = 1

for i in range(n_approx + 1):
    for j in range(n_approx + 1):
        # print('i = ', i)
        # print('j = ', j)
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
        # print('Часть от перемещения, мм', W_middle * 1000)
        equation_position += 1

Wn_lambda = lambdify([x, y], Wn)
W_middle = Wn_lambda(
    a / 2,
    b / 2
)
print('Перемещения в центре, мм', W_middle * 1000)
# W_edge = Wn_lambda(
#     0,
#     b / 2
# )
# print('Перемещения на краю, мм', W_edge * 1000)
# W_edge = Wn_lambda(
#     a / 2,
#     0
# )
# print('Перемещения на краю, мм', W_edge * 1000)
# W_edge = Wn_lambda(
#     a,
#     b
# )
# print('Перемещения на краю, мм', W_edge * 1000)

t2 = time.time()
t = t2 - t1
print('time', t)

# plot_graph(Wn_lambda, a, b)

