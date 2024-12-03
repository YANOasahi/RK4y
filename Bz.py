# import some libraries
import numpy as np

# import variables defined by myself
import variables_conditions as vc
import variables_Bz as vb

# import functions made by myself
import FourGaussian
import Enge
import position


# if magnet has trim coils, Bz : B0z * B-X curve * trim * fringe
# else if there is no trim coils, Bz : B0z * B-X curve * fringe
# else input parameter is wrong...
def Bz(trim, x, y):
    if trim == 'trim':
        Bz = vc.B0z *\
            (vb.para_twz[0]+vb.para_twz[1]*(x/1000) +
             vb.para_twz[2]*(x/1000)**2+vb.para_twz[3]*(x/1000)**3 +
             vb.para_twz[4]*(x/1000)**4+vb.para_twz[5]*(x/1000)**5 +
             vb.para_twz[6]*(x/1000)**6+vb.para_twz[7]*(x/1000)**7 +
             vb.para_twz[8]*(x/1000)**8)
        Bz = Bz+FourGaussian.Btrim_Sum(x)
        Bz = Bz*Enge.enge(y)*-1
        return Bz
    elif trim == 'no_trim':
        Bz = vc.B0z *\
            (vb.para_twz[0]+vb.para_twz[1]*(x/1000) +
             vb.para_twz[2]*(x/1000)**2+vb.para_twz[3]*(x/1000)**3 +
             vb.para_twz[4]*(x/1000)**4+vb.para_twz[5]*(x/1000)**5 +
             vb.para_twz[6]*(x/1000)**6+vb.para_twz[7]*(x/1000)**7 +
             vb.para_twz[8]*(x/1000)**8)
        Bz = Bz*Enge.enge(y)*-1
        return Bz
    else:
        print('You put a wrong input parameter!!')
        Bz = 0
        return Bz
