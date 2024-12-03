# import some libraries
import numpy as np
import matplotlib.pyplot as plt

# import variables defined by myself
import variables_FourGaussians as vf
import variables_tune as vt

# define gaussian


def Gaussian(x, amp, offset, mean, sigma):
    return amp*np.exp(-(x+offset-mean)**2/(2*sigma**2))

# sum up four gaussians
# after summing up, we have 10 arrays, Btrim


def Btrim_Sum(x):

    Btrim = np.zeros((10, len(x)))

    for i in range(10):
        for j in range(4):
            Btrim[i] += Gaussian(x, vf.amp[i][j], vf.offset[i][j],
                                 vf.mean[i][j], vf.sigma[i][j])

    # sum up 10 Btrims
    # after summing up, we have the magnetic field distribution of trim coils, Btrim_sum
    Btrim_sum = np.zeros(len(x))

    for i in range(10):
        Btrim_sum += Btrim[i]*vt.trim_current[i]/200

    return Btrim_sum

# check the shape of the function
x=np.arange(-300,300,0.5)
plt.plot(x,Btrim_Sum(x))
plt.show()
