import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos, solve, expand

from app.Functions import wn, eq1_integral, right_part_eq1, eq2_integral, right_part_eq2, eq3_integral, right_part_eq3, \
    eq4_integral, right_part_eq4

a = 5  # m
b = 6  # m
N1 = 0
E = 27500  # MPa
mu = 0.2
h = 0.2  # m
q = 0.01  # MPa

D = E * h ** 3 / (12 * (1 - mu ** 2))

t1 = time.time()

# Сбор основной системы уравнений
number_variables = 4  # wi
left_part: np.ndarray = np.zeros((number_variables, number_variables))
right_part: np.ndarray = np.zeros(number_variables)

i = 0
j = 0
w1_in_Wn, w2_in_Wn, w3_in_Wn, w4_in_Wn = wn(i, j, a, b)

left_part[0][0] = eq1_integral(i, j, a, b, D, q, w1_in_Wn)[0]
left_part[0][1] = eq1_integral(i, j, a, b, D, q, w2_in_Wn)[0]
left_part[0][2] = eq1_integral(i, j, a, b, D, q, w3_in_Wn)[0]
left_part[0][3] = eq1_integral(i, j, a, b, D, q, w4_in_Wn)[0]
right_part[0] = right_part_eq1(i, j, a, b, q)[0]
left_part[1][0] = eq2_integral(i, j, a, b, D, q, w1_in_Wn)[0]
left_part[1][1] = eq2_integral(i, j, a, b, D, q, w2_in_Wn)[0]
left_part[1][2] = eq2_integral(i, j, a, b, D, q, w3_in_Wn)[0]
left_part[1][3] = eq2_integral(i, j, a, b, D, q, w4_in_Wn)[0]
right_part[1] = right_part_eq2(i, j, a, b, q)[0]
left_part[2][0] = eq3_integral(i, j, a, b, D, q, w1_in_Wn)[0]
left_part[2][1] = eq3_integral(i, j, a, b, D, q, w2_in_Wn)[0]
left_part[2][2] = eq3_integral(i, j, a, b, D, q, w3_in_Wn)[0]
left_part[2][3] = eq3_integral(i, j, a, b, D, q, w4_in_Wn)[0]
right_part[2] = right_part_eq3(i, j, a, b, q)[0]
left_part[3][0] = eq4_integral(i, j, a, b, D, q, w1_in_Wn)[0]
left_part[3][1] = eq4_integral(i, j, a, b, D, q, w2_in_Wn)[0]
left_part[3][2] = eq4_integral(i, j, a, b, D, q, w3_in_Wn)[0]
left_part[3][3] = eq4_integral(i, j, a, b, D, q, w4_in_Wn)[0]
right_part[3] = right_part_eq4(i, j, a, b, q)[0]
print(left_part)
print(right_part)

vector_w: np.ndarray = np.linalg.solve(left_part, right_part)

x, y = symbols('x y')
w1, w2, w3, w4 = symbols('w1 w2 w3 w4')
Wn = (w1_in_Wn) * w1 + (w2_in_Wn) * w2 + (w3_in_Wn) * w3 + (w4_in_Wn) * w4
Wn_lambda = lambdify([w1, w2, w3, w4, x, y], Wn)
W_middle = Wn_lambda(
    vector_w[0],
    vector_w[1],
    vector_w[2],
    vector_w[3],
    a / 2,
    b / 2
)
print('Перемещения, мм', W_middle * 1000)

t2 = time.time()
t = t2 - t1
print('time', t)

