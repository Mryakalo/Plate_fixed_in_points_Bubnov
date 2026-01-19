import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos, solve

a = 5  # m
b = 6  # m
N1 = 0
E = 27500  # MPa
mu = 0.2
h = 0.2  # m
q = 0.01  # MPa

D = E * h ** 3 / (12 * (1 - mu ** 2))

# Сбор основной системы уравнений
number_variables = 4  # wi
left_part: np.ndarray = np.zeros((number_variables, number_variables))
right_part: np.ndarray = np.zeros(number_variables)

x, y = symbols('x y')
w1_00, w2_00, w3_00, w4_00 = symbols('w1_00 w2_00 w3_00 w4_00')
w1_01, w2_01, w3_01, w4_01 = symbols('w1_01 w2_01 w3_01 w4_01')
w1_11, w2_11, w3_11, w4_11 = symbols('w1_11 w2_11 w3_11 w4_11')
w1_10, w2_10, w3_10, w4_10 = symbols('w1_10 w2_10 w3_10 w4_10')
omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
# omega_lambda = lambdify([x, y], omega)
# print(omega)
# print(omega_lambda(a / 2, b / 2))

i = 0
j = 0
summa_00 = (w1_00 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w2_00 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b) +
            w3_00 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w4_00 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b))

i = 1
j = 0
summa_10 = (w1_10 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w2_10 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b) +
            w3_10 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w4_10 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b))

i = 1
j = 1
summa_11 = (w1_11 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w2_11 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b) +
            w3_11 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w4_11 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b))

i = 0
j = 1
summa_01 = (w1_01 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w2_01 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b) +
            w3_01 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
            w4_01 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b))

summa = simplify('(' + str(summa_00) + '+' + str(summa_10) + '+' + str(summa_11) + '+' + str(summa_01) + ')')

Wn_00 = simplify('(' + str(omega) + ') ** 2 * (' + str(summa_00) + ')')
Wn_10 = simplify('(' + str(omega) + ') ** 2 * (' + str(summa_10) + ')')
Wn_11 = simplify('(' + str(omega) + ') ** 2 * (' + str(summa_11) + ')')
Wn_01 = simplify('(' + str(omega) + ') ** 2 * (' + str(summa_01) + ')')

# print('summa_00 = ', summa_00)
# print('summa_10 = ', summa_10)
# print('summa_11 = ', summa_11)
# print('summa_01 = ', summa_01)

d4Wn_dx4 = '(' + str(diff(diff(diff(diff(Wn_00, x), x), x), x)) + ')'
d4Wn_dy4 = '(' + str(diff(diff(diff(diff(Wn_00, y), y), y), y)) + ')'
d4Wn_dx2_dy2 = '(' + str(diff(diff(diff(diff(Wn_00, x), x), y), y)) + ')'
nabla_Wn = simplify('(' + d4Wn_dx4 + '+' + d4Wn_dy4 + ' + 2 * ' + d4Wn_dx2_dy2 + ')')

i = 0
j = 0

eq1_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq1_00 = integrate(integrate(eq1_integrand.expand(), (x, 0, a)), (y, 0, b))
left_part[0][0] = diff(eq1_00, w1_00)
left_part[0][1] = diff(eq1_00, w2_00)
left_part[0][2] = diff(eq1_00, w3_00)
left_part[0][3] = diff(eq1_00, w4_00)
lamda_eq1 = lambdify([w1_00, w2_00, w3_00, w4_00], eq1_00)
right_part[0] = - lamda_eq1(0, 0, 0, 0)
print(left_part)

eq2_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq2_00 = integrate(integrate(eq2_integrand.expand(), (x, 0, a)), (y, 0, b))

eq3_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq3_00 = integrate(integrate(eq3_integrand.expand(), (x, 0, a)), (y, 0, b))

eq4_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq4_00 = integrate(integrate(eq4_integrand.expand(), (x, 0, a)), (y, 0, b))

i = 1
j = 0

eq1_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq1_10 = integrate(integrate(eq1_integrand.expand(), (x, 0, a)), (y, 0, b))

eq2_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq2_10 = integrate(integrate(eq2_integrand.expand(), (x, 0, a)), (y, 0, b))

eq3_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq3_10 = integrate(integrate(eq3_integrand.expand(), (x, 0, a)), (y, 0, b))

eq4_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq4_10 = integrate(integrate(eq4_integrand.expand(), (x, 0, a)), (y, 0, b))

i = 1
j = 1

eq1_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq1_11 = integrate(integrate(eq1_integrand.expand(), (x, 0, a)), (y, 0, b))

eq2_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq2_11 = integrate(integrate(eq2_integrand.expand(), (x, 0, a)), (y, 0, b))

eq3_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq3_11 = integrate(integrate(eq3_integrand.expand(), (x, 0, a)), (y, 0, b))

eq4_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq4_11 = integrate(integrate(eq4_integrand.expand(), (x, 0, a)), (y, 0, b))

i = 0
j = 1

eq1_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq1_01 = integrate(integrate(eq1_integrand.expand(), (x, 0, a)), (y, 0, b))

eq2_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq2_01 = integrate(integrate(eq2_integrand.expand(), (x, 0, a)), (y, 0, b))

eq3_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq3_01 = integrate(integrate(eq3_integrand.expand(), (x, 0, a)), (y, 0, b))

eq4_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq4_01 = integrate(integrate(eq4_integrand.expand(), (x, 0, a)), (y, 0, b))

system = [eq1_00, eq2_00, eq3_00, eq4_00,
          eq1_10, eq2_10, eq3_10, eq4_10,
          eq1_11, eq2_11, eq3_11, eq4_11,
          eq1_01, eq2_01, eq3_01, eq4_01]

solutions = solve(system, w1_00, w2_00, w3_00, w4_00,
                  w1_10, w2_10, w3_10, w4_10,
                  w1_11, w2_11, w3_11, w4_11,
                  w1_01, w2_01, w3_01, w4_01,dict=True)
# print(solutions[0][w1])

Wn_lambda = lambdify([w1_00, w2_00, w3_00, w4_00,
                           w1_10, w2_10, w3_10, w4_10,
                           w1_11, w2_11, w3_11, w4_11,
                           w1_01, w2_01, w3_01, w4_01, x, y], Wn)
W_middle = Wn_lambda(
    solutions[0][w1_00],
    solutions[0][w2_00],
    solutions[0][w3_00],
    solutions[0][w4_00],
    solutions[0][w1_10],
    solutions[0][w2_10],
    solutions[0][w3_10],
    solutions[0][w4_10],
    solutions[0][w1_11],
    solutions[0][w2_11],
    solutions[0][w3_11],
    solutions[0][w4_11],
    solutions[0][w1_01],
    solutions[0][w2_01],
    solutions[0][w3_01],
    solutions[0][w4_01],
    a / 2,
    b / 2
)
print('Перемещения, мм', W_middle * 1000)