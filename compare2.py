import numpy as np
import matplotlib.pyplot as plt

plt.style.use('./ticksin_helvetica.mplstyle')

box1 = plt.figure(figsize=(10, 6))
fig1_1 = box1.add_subplot(1, 1, 1)

data=np.genfromtxt('./bzz_yano.dat')
fig1_1 = plt.plot(data[:,0], data[:,1], label="Yano results")

data2=np.genfromtxt('./bzz_abesan.dat')
fig1_1 = plt.plot(data2[:,0], data2[:,1], label="Abe-san results")

plt.legend()
plt.show()