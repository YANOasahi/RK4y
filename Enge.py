# importing some libraries
import numpy as np
import matplotlib.pyplot as plt

# import variables made bu myself
import variables_Enge

# define 6order Enge function
def Enge(x):
    return 1/(1+(np.exp(variables_Enge.enge[0]+(x*variables_Enge.enge[1])+\
                (x**2*variables_Enge.enge[2])+(x**3*variables_Enge.enge[3])+\
                (x**4*variables_Enge.enge[4])+(x**5*variables_Enge.enge[5]))))



# check the shape of the function
x = np.arange(-10, 10, 0.5)

plt.plot(x, Enge(x))
plt.show()