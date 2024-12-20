import numpy as np
import matplotlib.pyplot as plt

# データの読み込み
data = np.genfromtxt('./by_map.dat')

x_axis = np.array(data[:, 0])
y_axis = np.array(data[:, 1])
z_axis = np.array(data[:, 2])
by = np.array(data[:, 3])

# 3Dスキャッタープロット
box1 = plt.figure(figsize=(14, 6))

# x-y 平面での3Dスキャッター
fig1_1 = box1.add_subplot(1, 3, 1, projection='3d')
fig1_1.scatter(x_axis, y_axis, by, c=by, cmap='viridis', s=10)
fig1_1.set_xlabel('X-axis')
fig1_1.set_ylabel('Y-axis')
fig1_1.set_zlabel('BY')
fig1_1.set_title('X-Y vs BY')

# y-z 平面での3Dスキャッター
fig1_2 = box1.add_subplot(1, 3, 2, projection='3d')
fig1_2.scatter(y_axis, z_axis, by, c=by, cmap='viridis', s=10)
fig1_2.set_xlabel('Y-axis')
fig1_2.set_ylabel('Z-axis')
fig1_2.set_zlabel('BY')
fig1_2.set_title('Y-Z vs BY')

# x-z 平面での3Dスキャッター
fig1_3 = box1.add_subplot(1, 3, 3, projection='3d')
fig1_3.scatter(x_axis, z_axis, by, c=by, cmap='viridis', s=10)
fig1_3.set_xlabel('X-axis')
fig1_3.set_ylabel('Z-axis')
fig1_3.set_zlabel('BY')
fig1_3.set_title('X-Z vs BY')

plt.tight_layout()
plt.show()