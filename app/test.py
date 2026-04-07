from scipy import integrate
from sympy import lambdify, symbols, diff, simplify, sin, cos, exp, expand

x, y = symbols('x y')
a = 6  # m
b = 6  # m
integrand = 1.43537272989276*(sin(0.523598775598299*x))**3*sin(0.523598775598299*y) + 2.87074545978553*(sin(0.523598775598299*x))**2*(sin(0.523598775598299*y))**2 + 1.43537272989276*sin(0.523598775598299*x)*(sin(0.523598775598299*y))**3
integrand_lambda = lambdify([x, y], integrand)
coef_wi_in_equation = integrate.nquad(integrand_lambda, [[0, a], [0, b]])
print('integrand = ', integrand, '   coef_wi_in_equation = ', coef_wi_in_equation)
