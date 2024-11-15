# importing some libraries
import numpy as np
import matplotlib.pyplot as plt

# define Enge function
def Enge(x,first,second,third,fourth):
    return 1/(1+(np.exp(first+(x*second)+(x**2*third)+(x**3*fourth))))

# check the shape of the function
x=np.arange(-10,10,0.5)

plt.plot(x,Enge(x,0.95,1,0.2,0.05))
plt.show()