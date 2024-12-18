import numpy as np
import matplotlib.pyplot as plt

box1 = plt.figure(figsize=(7, 6))
fig1_1 = box1.add_subplot(1, 1, 1)

yano=np.genfromtxt('./rk45_output_max1ps.dat')
# yano=np.genfromtxt('./rk45_output.dat')
fig1_1 = plt.plot(yano[:,0], yano[:,1], label="Yano results")
abesan=np.genfromtxt('./kidou_long.dat')
fig1_1 = plt.plot(abesan[:,1], abesan[:,2], label="Abe-san results")

plt.legend()
plt.show()