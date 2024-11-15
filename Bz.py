# import some libraries
import numpy as np

# import variables defined by myself
import variables_tune
import variables_Bz

# import functions made by myself
import FourGaussian
import Enge


# if magnet has trim coils, Bz : B0z * B-X curve * trim * fringe
# else if there is no trim coils, Bz : B0z * B-X curve * fringe
# else input parameter is wrong...
def Bx(trim, x):
    if trim == 'trim':
        Bz = variables_tune.B0z *\
            (variables_Bz.para_wz[0]+variables_Bz.para_wz[1]*x +
             variables_Bz.para_wz[2]*x**2+variables_Bz.para_wz[3]*x**3 +
             variables_Bz.para_wz[4]*x**4+variables_Bz.para_wz[5]*x**5 +
             variables_Bz.para_wz[6]*x**6+variables_Bz.para_wz[7]*x**7 +
             variables_Bz.para_wz[8]*x**8+variables_Bz.para_wz[9]*x**9) +\
            FourGaussian.Btrim_Sum(x) *\
            Enge.enge(x)*-1
        return Bz
    elif trim == 'no_trim':
        Bz = variables_tune.B0z *\
            (variables_Bz.para_wz[0]+variables_Bz.para_wz[1]*x +
             variables_Bz.para_wz[2]*x**2+variables_Bz.para_wz[3]*x**3 +
             variables_Bz.para_wz[4]*x**4+variables_Bz.para_wz[5]*x**5 +
             variables_Bz.para_wz[6]*x**6+variables_Bz.para_wz[7]*x**7 +
             variables_Bz.para_wz[8]*x**8+variables_Bz.para_wz[9]*x**9) *\
            Enge.enge(x)*-1
        return Bz
    else:
        print('You put a wrong input parameter!!')
        Bz = 0
        return Bz
