from sympy import symbols, diff

def calc_moment(
        D,
        mu,
        Wn,

):
    x, y, z = symbols('x y z')

    dWn_dx = diff(Wn, x)
    d2Wn_dx = diff(dWn_dx, x)
    dWn_dy = diff(Wn, y)
    d2Wn_dy = diff(dWn_dy, y)

    # print('d2Wn_dx = ', d2Wn_dx)
    # print('d2Wn_dy = ', d2Wn_dy)

    # sigma_x = - E * z / (1 - mu) * (d2Wn_dx + mu * (d2Wn_dy))
    # sigma_y = - E * z / (1 - mu) * (d2Wn_dy + mu * (d2Wn_dx))

    Mx = str(- D) + ' * (' + str(d2Wn_dx) + ' + ' + str(mu) + ' * (' + str(d2Wn_dy) + '))'
    My = str(- D) + ' * (' + str(d2Wn_dy) + ' + ' + str(mu) + ' * (' + str(d2Wn_dx) + '))'
    # print('Mx = ', Mx)

    return Mx, My