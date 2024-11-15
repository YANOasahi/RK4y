# import some libraries
import numpy as np
import matplotlib.pyplot as plt

# import variables defined by myself
import variables

# define gaussian
def Gaussian(x, amp, offset, mean, sigma):
    return amp*np.exp(-(x+offset-mean)**2/(2*sigma**2))

# sum up four gaussians
# after summing up, we have 10 arrays, Btrim
def Btrim_Sum(x):
    
    Btrim = np.zeros(10)

    for i in range(10):
        for j in range(4):
            Btrim[i] += Gaussian(x, variables.amp[i][j], variables.offset[i][j], 
                                 variables.mean[i][j], variables.sigma[i][j])

    # sum up 10 Btrims
    # after summing up, we have the magnetic field distribution of trim coils, Btrim_sum
    Btrim_sum=np.zeros()

    for i in range(10):
        Btrim_sum+=Btrim[i]

    # # check the shape of the function
    # for i in range(10):
    #      plt.plot(x,Btrim[i])

    # plt.plot(x,Btrim_sum+0.05)

    # plt.show()

    return Btrim_sum