import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, simplify, sin, cos, expand
from scipy import integrate

def omega_polynomial(a, b):
    x, y = symbols('x y')
    return simplify('(' + str(sin(m.pi * x / a) * (sin(m.pi * y / b) - 2) + sin(m.pi * y / b) * (sin(m.pi * x / a) - 2)) + ')')

def wn_polynomial(i, j, a, b):
    x, y = symbols('x y')

    coef_w1_in_sum = cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)

    coef_w1_in_wn = simplify('(' + str(omega_polynomial(a, b)) + ') * (' + str(coef_w1_in_sum) + ')')

    return coef_w1_in_wn


def nabla_wn_polynomial(w_in_wn):
    x, y = symbols('x y')
    d4wn_dx4 = '(' + str(diff(expand(diff(expand(diff(expand(diff(expand(w_in_wn), x)), x)), x)), x)) + ')'
    d4wn_dy4 = '(' + str(diff(expand(diff(expand(diff(expand(diff(expand(w_in_wn), y)), y)), y)), y)) + ')'
    d4wn_dx2_dy2 = '(' + str(diff(expand(diff(expand(diff(expand(diff(expand(w_in_wn), x)), x)), y)), y)) + ')'
    nabla = simplify('(' + d4wn_dx4 + '+' + d4wn_dy4 + ' + 2 * (' + d4wn_dx2_dy2 + '))')

    return nabla


def eq1_integral_polynomial(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn_polynomial(w_in_wn)) * omega_polynomial(a, b) * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return coef_wi_in_equation


def right_part_eq1_polynomial(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega_polynomial(a, b) * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation



