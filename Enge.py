# importing some libraries
import numpy as np
import matplotlib.pyplot as plt
from math import fabs

# import variables made bu myself
import variables_Enge as ve


# define 6order Enge function
def enge(input):
    length=525
    aperture=80
    x=(fabs(input)-length)/aperture
    expo = ve.enge[0]+x*ve.enge[1]+(x**2)*ve.enge[2]+ \
           (x**3)*ve.enge[3]+(x**4)*ve.enge[4]+(x**5)*ve.enge[5]
    denominator = 1+np.exp(expo)
    return 1 / denominator


# # check the shape of the function
# x = np.arange(-1000, 1000, 0.1)
# plt.plot(x, enge(x))
# plt.show()
