import numpy as np
import math as m
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos, solve

a = 10  # m
b = 20  # m
N1 = 0
E = 30000  # MPa
mu = 0.2
h = 0.2  # m
q = 0.01  # MPa

D = E * h ** 3 / (12 * (1 - mu ** 2))

x, y = symbols('x y')
w1, w2, w3, w4 = symbols('w1 w2 w3 w4')
omega = sin(m.pi * x / a) + sin(m.pi * y / b)

i = 0
j = 0
summa = (w1 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
         w2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b) +
         w3 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b) +
         w4 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b))

Wn = omega ** 2 * summa

nabla_W = diff(diff(Wn, x), y)

eq1_integrand = (D * nabla_W - q) * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq1 = integrate(integrate(eq1_integrand, (x, 0, a)), (y, 0, b))

eq2_integrand = (D * nabla_W - q) * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq2 = integrate(integrate(eq2_integrand, (x, 0, a)), (y, 0, b))

eq3_integrand = (D * nabla_W - q) * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
eq3 = integrate(integrate(eq3_integrand, (x, 0, a)), (y, 0, b))

eq4_integrand = (D * nabla_W - q) * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
eq4 = integrate(integrate(eq4_integrand, (x, 0, a)), (y, 0, b))

system = [eq1, eq2, eq3, eq4]
solutions = solve(system, w1, w2, w3, w4)
print(solutions[0])

Wn_lambda = lambdify([w1, w2, w3, w4, x, y], Wn)
# W_middle = Wn_lambda(solutions[0])