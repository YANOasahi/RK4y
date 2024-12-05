# importing some libraries
import numpy as np
import matplotlib.pyplot as plt

# import variables made bu myself
import variables_Enge as ve

# define 6order Enge function
def enge(x):
    if x<10:
        return 1/(1+(np.exp(ve.enge[0]+(x*ve.enge[1])+\
                ((x**2)*ve.enge[2])+((x**3)*ve.enge[3])+\
                ((x**4)*ve.enge[4])+((x**5)*ve.enge[5]))))
    else:
        return 0.0


# def enge(x):
#     return 1/(1+(np.exp(ve.enge[0]+(x*ve.enge[1]) +
#                         ((x**2)*ve.enge[2])+((x**3)*ve.enge[3]) +
#                         ((x**4)*ve.enge[4])+((x**5)*ve.enge[5]))))


# # check the shape of the function
# x = np.arange(-10, 10, 0.1)
# plt.plot(x, enge(x))
# plt.show()
