import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos, solve, expand

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
w1, w2, w3, w4 = symbols('w1 w2 w3 w4')
omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
# omega_lambda = lambdify([x, y], omega)
# print(omega)
# print(omega_lambda(a / 2, b / 2))

i = 0
j = 0
summa = (w1 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
         w2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b) +
         w3 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
         w4 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b))

Wn = simplify('(' + str(omega) + ') ** 2 * (' + str(summa) + ')')

d4Wn_dx4 = '(' + str(diff(diff(diff(diff(Wn, x), x), x), x)) + ')'
d4Wn_dy4 = '(' + str(diff(diff(diff(diff(Wn, y), y), y), y)) + ')'
d4Wn_dx2_dy2 = '(' + str(diff(diff(diff(diff(Wn, x), x), y), y)) + ')'
# print(d4Wn_dx2_dy2)
nabla_Wn = simplify('(' + d4Wn_dx4 + '+' + d4Wn_dy4 + ' + 2 * ' + d4Wn_dx2_dy2 + ')')

print('nabla_Wn =', nabla_Wn)
print('D * nabla_Wn =', D * nabla_Wn)
print('omega ** 2 =', omega ** 2)

eq1_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
print('eq1_integrand =', eq1_integrand.expand())

print('hhh', diff(eq1_integrand, w1))

eq1_integrand_expand = eq1_integrand.expand()
t1 = time.time()
coef_w1 = diff(eq1_integrand_expand, w1)
left_part[0][0] = integrate(integrate(coef_w1, (x, 0, a)), (y, 0, b))
coef_w2 = diff(eq1_integrand_expand, w2)
left_part[0][1] = integrate(integrate(coef_w2, (x, 0, a)), (y, 0, b))
coef_w3 = diff(eq1_integrand_expand, w3)
left_part[0][2] = integrate(integrate(coef_w3, (x, 0, a)), (y, 0, b))
coef_w4 = diff(eq1_integrand_expand, w4)
left_part[0][3] = integrate(integrate(coef_w4, (x, 0, a)), (y, 0, b))
free_member = - (eq1_integrand_expand - coef_w1 * w1 - coef_w2 * w2 - coef_w3 * w3 - coef_w4 * w4)
print(integrate(integrate(free_member, (x, 0, a)), (y, 0, b)))
# right_part[0] = integrate(integrate(free_member, (x, 0, a)), (y, 0, b))
t2 = time.time()
t = t2 - t1
print('time', t)
print(left_part)
print(right_part)
t1 = time.time()
eq1 = integrate(integrate(eq1_integrand_expand, (x, 0, a)), (y, 0, b))
left_part[0][0] = diff(eq1, w1)
left_part[0][1] = diff(eq1, w2)
left_part[0][2] = diff(eq1, w3)
left_part[0][3] = diff(eq1, w4)
right_part[0] = - (eq1 - left_part[0][0] * w1 - left_part[0][1] * w2 - left_part[0][2] * w3 - left_part[0][3] * w4)
t2 = time.time()
t = t2 - t1
print('time', t)
print(left_part)
print(right_part)

print('eq1 =', eq1)

eq2_integrand = (D * nabla_Wn - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
print('eq2_integrand =', eq2_integrand.expand())
eq2 = integrate(integrate(eq2_integrand.expand(), (x, 0, a)), (y, 0, b))
print('eq2 =', eq2)

eq3_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
print('eq3_integrand =', eq3_integrand.expand())
eq3 = integrate(integrate(eq3_integrand.expand(), (x, 0, a)), (y, 0, b))
print('eq3 =', eq3)

eq4_integrand = (D * nabla_Wn - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
print('eq4_integrand =', eq4_integrand.expand())
eq4 = integrate(integrate(eq4_integrand.expand(), (x, 0, a)), (y, 0, b))
print('eq4 =', eq4)

system = [eq1, eq2, eq3, eq4]
solutions = solve(system, w1, w2, w3, w4, dict=True)
print(solutions[0])

Wn_lambda = lambdify([w1, w2, w3, w4, x, y], Wn)
W_middle = Wn_lambda(
    solutions[0][w1],
    solutions[0][w2],
    solutions[0][w3],
    solutions[0][w4],
    a / 2,
    b / 2
)
print('Перемещения, мм', W_middle * 1000)

