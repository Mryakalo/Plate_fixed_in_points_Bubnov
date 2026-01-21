import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos
from scipy import integrate

def omega(a, b):
    x, y = symbols('x y')
    return simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b)) + ')')
    # return simplify('(' + str((a - x) ** 2 * x ** 2 + (b - y) ** 2 * y ** 2) + ')')

def wn(i, j, a, b):
    x, y = symbols('x y')

    coef_w1_in_sum = sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    coef_w2_in_sum = sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    coef_w3_in_sum = cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    coef_w4_in_sum = cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)

    coef_w1_in_wn = simplify('(' + str(omega(a, b)) + ') ** 2 * (' + str(coef_w1_in_sum) + ')')
    coef_w2_in_wn = simplify('(' + str(omega(a, b)) + ') ** 2 * (' + str(coef_w2_in_sum) + ')')
    coef_w3_in_wn = simplify('(' + str(omega(a, b)) + ') ** 2 * (' + str(coef_w3_in_sum) + ')')
    coef_w4_in_wn = simplify('(' + str(omega(a, b)) + ') ** 2 * (' + str(coef_w4_in_sum) + ')')

    return coef_w1_in_wn, coef_w2_in_wn, coef_w3_in_wn, coef_w4_in_wn


def nabla_wn(w_in_wn):
    x, y = symbols('x y')
    d4wn_dx4 = '(' + str(diff(diff(diff(diff(w_in_wn, x), x), x), x)) + ')'
    d4wn_dy4 = '(' + str(diff(diff(diff(diff(w_in_wn, y), y), y), y)) + ')'
    d4wn_dx2_dy2 = '(' + str(diff(diff(diff(diff(w_in_wn, x), x), y), y)) + ')'
    nabla = simplify('(' + d4wn_dx4 + '+' + d4wn_dy4 + ' + 2 * (' + d4wn_dx2_dy2 + '))')

    return nabla


def eq1_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return coef_wi_in_equation


def eq2_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return coef_wi_in_equation


def eq3_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return coef_wi_in_equation


def eq4_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return coef_wi_in_equation


def right_part_eq1(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) ** 2 * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation


def right_part_eq2(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) ** 2 * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation


def right_part_eq3(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) ** 2 * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation


def right_part_eq4(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) ** 2 * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation
