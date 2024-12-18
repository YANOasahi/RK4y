# importing some libraries
import numpy as np
import matplotlib.pyplot as plt
from math import fabs
from fractions import Fraction

# import variables made bu myself
import variables as var


# define 6order Enge function
def enge(input):
    length = 525.0
    aperture = 80.0
    x = (fabs(input)-length)/aperture
    expo = Fraction(var.enge[0]+x*var.enge[1]+np.square(x)*var.enge[2] +
                    np.power(x,3)*var.enge[3]+np.power(x,4)*var.enge[4]+np.power(x,5)*var.enge[5])
    denominator = Fraction(1+np.exp(float(expo)))
    return np.reciprocal(float(denominator))  # returun 1/denominator


# # check the shape of the function
# x = np.arange(-1000, 1000, 0.1)
# plt.plot(x, enge(x))
# plt.show()
