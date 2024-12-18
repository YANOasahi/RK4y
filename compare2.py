import numpy as np
import matplotlib.pyplot as plt

plt.style.use('./ticksin_helvetica.mplstyle')

box1 = plt.figure(figsize=(10, 6))
fig1_1 = box1.add_subplot(1, 1, 1)

# data1=np.genfromtxt('./bzz_yano.dat')
# fig1_1 = plt.plot(data1[:,0], data1[:,1], label="Yano 1e-6")

data2=np.genfromtxt('./bzz_yano2.dat')
fig1_1 = plt.plot(data2[:,0]*1000, data2[:,1], label="Yano")

# remove dupulication
input_file = "bzz_abesan.dat"  # origin
output_file = "bzz_abesan2.dat"  # new output

with open(input_file, "r") as infile:
    unique_lines = sorted(set(line.strip() for line in infile))

with open(output_file, "w") as outfile:
    outfile.write("\n".join(unique_lines) + "\n")

data3=np.genfromtxt('./bzz_abesan2.dat')
fig1_1 = plt.plot(data3[:,0], data3[:,1], label="Abe-san")

plt.legend()
plt.show()