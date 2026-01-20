import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos, solve, expand

def wn(i, j, a, b):
    x, y = symbols('x y')
    omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')

    coef_w1_in_sum = sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    coef_w2_in_sum = sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    coef_w3_in_sum = cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    coef_w4_in_sum = cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)

    coef_w1_in_wn = simplify('(' + str(omega) + ') ** 2 * (' + str(coef_w1_in_sum) + ')')
    coef_w2_in_wn = simplify('(' + str(omega) + ') ** 2 * (' + str(coef_w2_in_sum) + ')')
    coef_w3_in_wn = simplify('(' + str(omega) + ') ** 2 * (' + str(coef_w3_in_sum) + ')')
    coef_w4_in_wn = simplify('(' + str(omega) + ') ** 2 * (' + str(coef_w4_in_sum) + ')')

    return coef_w1_in_wn, coef_w2_in_wn, coef_w3_in_wn, coef_w4_in_wn


def nabla_wn(w_in_wn):

    x, y = symbols('x y')
    d4wn_dx4 = '(' + str(diff(diff(diff(diff(w_in_wn, x), x), x), x)) + ')'
    d4wn_dy4 = '(' + str(diff(diff(diff(diff(w_in_wn, y), y), y), y)) + ')'
    d4wn_dx2_dy2 = '(' + str(diff(diff(diff(diff(w_in_wn, x), x), y), y)) + ')'
    nabla = simplify('(' + d4wn_dx4 + '+' + d4wn_dy4 + ' + 2 * ' + d4wn_dx2_dy2 + ')')

    return nabla


def eq1_integral(i, j, a, b, D, q, w_in_wn):

    x, y = symbols('x y')
    omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
    integrand = (D * nabla_wn(w_in_wn) - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    eq1_integrand_expand = integrand.expand()
    eq1 = integrate(integrate(eq1_integrand_expand, (x, 0, a)), (y, 0, b))

    return eq1


def eq2_integral(i, j, a, b, D, q, w_in_wn):

    x, y = symbols('x y')
    omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
    integrand = (D * nabla_wn(w_in_wn) - q) * omega ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    eq2_integrand_expand = integrand.expand()
    eq2 = integrate(integrate(eq2_integrand_expand, (x, 0, a)), (y, 0, b))

    return eq2


def eq3_integral(i, j, a, b, D, q, w_in_wn):

    x, y = symbols('x y')
    omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
    integrand = (D * nabla_wn(w_in_wn) - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    eq3_integrand_expand = integrand.expand()
    eq3 = integrate(integrate(eq3_integrand_expand, (x, 0, a)), (y, 0, b))

    return eq3


def eq4_integral(i, j, a, b, D, q, w_in_wn):

    x, y = symbols('x y')
    omega = simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
    integrand = (D * nabla_wn(w_in_wn) - q) * omega ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    eq4_integrand_expand = integrand.expand()
    eq4 = integrate(integrate(eq4_integrand_expand, (x, 0, a)), (y, 0, b))

    return eq4

