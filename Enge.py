# importing some libraries
import numpy as np
import matplotlib.pyplot as plt
from math import fabs
from fractions import Fraction

# import variables made bu myself
import variables_Enge as ve


# define 6order Enge function
def enge(input):
    length=525.0
    aperture=80.0
    x=(fabs(input)-length)/aperture
    expo = Fraction(ve.enge[0]+x*ve.enge[1]+(x**2)*ve.enge[2]+ \
           (x**3)*ve.enge[3]+(x**4)*ve.enge[4]+(x**5)*ve.enge[5])
    denominator = Fraction(1+np.exp(float(expo)))
    return 1 / denominator


# # check the shape of the function
# x = np.arange(-1000, 1000, 0.1)
# plt.plot(x, enge(x))
# plt.show()
