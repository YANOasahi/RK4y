import numpy as np

import variables_conditions as vc

import Bz


def diff_px(x, y, v, gamma, trim):
    return x*v*gamma*Bz.bz(trim, x, y)*np.sqrt(x**2+y**2)


def diff_py(x, y, v, gamma, trim):
    return y*v*gamma*Bz.bz(trim, x, y)*np.sqrt(x**2+y**2)


def diff_x(px):
    return px/vc.mass

def diff_y(py):
    return py/vc.mass