# import some libraries
import numpy as np
import matplotlib.pyplot as plt

# import variables defined by myself
import variables_FourGaussians as vf
import variables_conditions as vc

# define gaussian


def Gaussian(x, amp, offset, mean, sigma):
    return amp*np.exp(-(x+offset-mean)**2/(2*sigma**2))

# sum up four gaussians
# after summing up, we have 10 arrays, Btrim


def Btrim_Sum(x):
    Btrim = np.zeros(10)
    for i in range(10):
        Btrim[i] = sum(Gaussian(x, vf.amp[i][j], vf.offset[i][j],
                       vf.mean[i][j], vf.sigma[i][j]) for j in range(4))*vc.trim_current[i]/200
    Btrim_sum = np.sum(Btrim)
    return Btrim_sum


# # check the shape of the function
# x = np.arange(-300, 300, 0.5)
# plt.plot(x, [Btrim_Sum(val) for val in x])
# plt.show()
