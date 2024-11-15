# importing some libraries
import numpy as np
import matplotlib.pyplot as plt

# define Enge function
def FourGaussian(x,amp,sigma,mean1,mean2,mean3,mean4):
    return amp*np.exp(-(x-mean1)**2/(sigma**2))+amp*np.exp(-(x-mean2)**2/(sigma**2))+amp*np.exp(-(x-mean3)**2/(sigma**2))+amp*np.exp(-(x-mean4)**2/(sigma**2))

# check the shape of the function
x=np.arange(-50,150,0.5)

plt.plot(x,FourGaussian(x,100,20,15,35,55,75))
plt.show()