import numpy as np
import math as m
import time
from sympy import lambdify, symbols, diff, integrate, simplify, sin, cos, exp, expand
from scipy import integrate

def omega(a, b):
    x, y = symbols('x y')
    return simplify('(' + str(sin(m.pi * x / a) + sin(m.pi * y / b) - 0.2 * sin(m.pi * y / b) * sin(m.pi * x / a)) + ')')
    # return simplify('(' + str((a - x) ** 2 * x ** 2 + (b - y) ** 2 * y ** 2) + ')')
    # return simplify('(' + str((a - x) * x + (b - y) * y) + ')')

def wn(i, j, a, b):
    x, y = symbols('x y')

    coef_w1_in_sum = sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    coef_w2_in_sum = sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    coef_w3_in_sum = cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    coef_w4_in_sum = cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)

    coef_w1_in_wn = simplify('(' + str(omega(a, b) * coef_w1_in_sum) + ')')
    coef_w2_in_wn = simplify('(' + str(omega(a, b) * coef_w2_in_sum) + ')')
    coef_w3_in_wn = simplify('(' + str(omega(a, b) * coef_w3_in_sum) + ')')
    coef_w4_in_wn = simplify('(' + str(omega(a, b) * coef_w4_in_sum) + ')')

    return coef_w1_in_wn, coef_w2_in_wn, coef_w3_in_wn, coef_w4_in_wn


def nabla_wn(w_in_wn):
    t1 = time.time()
    x, y = symbols('x y')
    d4wn_dx4 = '(' + str(diff(expand(diff(expand(diff(expand(diff(expand(w_in_wn), x)), x)), x)), x)) + ')'
    d4wn_dy4 = '(' + str(diff(expand(diff(expand(diff(expand(diff(expand(w_in_wn), y)), y)), y)), y)) + ')'
    d4wn_dx2_dy2 = '(' + str(diff(expand(diff(expand(diff(expand(diff(expand(w_in_wn), x)), x)), y)), y)) + ')'
    nabla = simplify('((' + d4wn_dx4 + ')+(' + d4wn_dy4 + ' )+ 2 * (' + d4wn_dx2_dy2 + '))')
    # print('w_in_wn = ', w_in_wn)
    # print('d4wn_dx4 = ', d4wn_dx4,
    #       ' d4wn_dy4 = ', d4wn_dy4,
    #       ' d4wn_dx2_dy2 = ', d4wn_dx2_dy2,
    #       ' nabla = ', nabla)
    t2 = time.time()
    t = t2 - t1
    # print('time diff = ', t)

    return nabla


def eq1_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = expand((D * nabla_wn(w_in_wn)) * omega(a, b) * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b))
    integrand_lambda = lambdify([x, y], integrand)
    # print('w_in_wn = ', w_in_wn)
    # print('nabla_wn(w_in_wn) = ', nabla_wn(w_in_wn))
    t1 = time.time()
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])
    # print('integrand = ', integrand, '   coef_wi_in_equation = ', coef_wi_in_equation)
    t2 = time.time()
    t = t2 - t1
    # print('time integr eq1 = ', t)

    return coef_wi_in_equation[0]


def eq2_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    t1 = time.time()
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])
    t2 = time.time()
    t = t2 - t1
    # print('time integr eq2 = ', t)

    return coef_wi_in_equation[0]


def eq3_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    t1 = time.time()
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])
    t2 = time.time()
    t = t2 - t1
    # print('time integr eq3 = ', t)

    return coef_wi_in_equation[0]


def eq4_integral(i, j, a, b, D, w_in_wn):
    x, y = symbols('x y')
    integrand = (D * nabla_wn(w_in_wn)) * omega(a, b) * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    t1 = time.time()
    coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])
    t2 = time.time()
    t = t2 - t1
    # print('time integr eq4 = ', t)

    return coef_wi_in_equation[0]


def right_part_eq1(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) * sin((2 * i + 1) * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation


def right_part_eq2(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) * sin((2 * i + 1) * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation


def right_part_eq3(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) * cos(2 * i * m.pi * x / a) * sin((2 * j + 1) * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation


def right_part_eq4(i, j, a, b, q):
    x, y = symbols('x y')
    integrand = q * omega(a, b) * cos(2 * i * m.pi * x / a) * cos(2 * j * m.pi * y / b)
    integrand_lambda = lambdify([x, y], integrand)
    free_member_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])

    return free_member_in_equation
