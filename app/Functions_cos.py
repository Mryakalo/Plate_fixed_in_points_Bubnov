import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos
from scipy import integrate

def omega_cos(a, b):
    x, y = symbols('x y')
    return simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')

def wn_cos(i, j, a, b):
    x, y = symbols('x y')

    coef_w1_in_sum = 3.14 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)

    coef_w1_in_wn = simplify('(' + str(omega_cos(a, b)) + ') ** 2 * (' + str(coef_w1_in_sum) + ')')

    return coef_w1_in_wn


def nabla_wn_cos(w_in_wn):
    x, y = symbols('x y')
    d4wn_dx4 = '(' + str(diff(diff(diff(diff(w_in_wn, x), x), x), x)) + ')'
    d4wn_dy4 = '(' + str(diff(diff(diff(diff(w_in_wn, y), y), y), y)) + ')'
    d4wn_dx2_dy2 = '(' + str(diff(diff(diff(diff(w_in_wn, x), x), y), y)) + ')'
    nabla = simplify('(' + d4wn_dx4 + '+' + d4wn_dy4 + ' + 2 * (' + d4wn_dx2_dy2 + '))')

    return nabla


def eq1_integral_cos(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn_cos(w_in_wn)) * omega_cos(a, b) ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return coef_wi_in_equation


def right_part_eq1_cos(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega_cos(a, b) ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation



