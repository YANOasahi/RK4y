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
    # buf is a variable for temporary use
    buf = 0.0
    # Btrim is a list of 10 trim coils
    Btrim = [(10)]
    # Btrim_sum is used in Bz.py
    Btrim_sum = 0.0

    for i in range(10):
        for j in range(4):
            buf += Gaussian(x, vf.amp[i][j], vf.offset[i][j],
                            vf.mean[i][j], vf.sigma[i][j])
            if j == 3:
                Btrim.append(buf)
                Btrim_sum += buf*vc.trim_current[i]/200

    return Btrim_sum


# check the shape of the function
# x = np.arange(-300, 300, 0.5)
# plt.plot(x, Btrim_Sum(x))
# plt.show()
