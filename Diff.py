import numpy as np

import variables_conditions as vc

import Bz


def diff_vx(x, y, vy, gamma):
    # print(f'factor of vb is {vc.z*vc.c/(vc.mass*vc.amu*gamma)}')
    return vc.z*vc.c*vy*Bz.BforXplane(x, y)/(vc.mass*vc.amu*gamma)


def diff_vy(x, y, vx, gamma):
    return vc.z*vc.c*-vx*Bz.BforXplane(x, y)/(vc.mass*vc.amu*gamma)


def diff_x(vx):
    return vx


def diff_y(vy):
    return vy
