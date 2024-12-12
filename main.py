import numpy as np
import matplotlib.pyplot as plt
import time
import variables_position as vp
import variables_conditions as vc
import Diff
import Bz
import position

# for measuring execution time
start = time.perf_counter()

#  set brho
j = 0.0
dp = 0.1 * j
brho = vc.brho0 * (1 + dp / 100)

# set emittance
l = 0.0
emi = l * 20  # 20 steps
a_max = np.sqrt(emi / vc.betax) * 0.95
a_bin = 2 * a_max / 15
k = 0.0
a_init = k * a_bin - a_max

# set dispersion [m/%] in the ring
dispersion = -3.3494 * dp**6 + 2.2088 * dp**5 + 1.5699 * dp**4 - \
    0.731 * dp**3 - 0.1573 * dp**2 + 0.0121 * dp + 7.0737

# q/m
MoQ = vc.mass / vc.z
# momentum based on brho
mom = vc.c * vc.z * brho * 1E-6
# beam energy based on momentum
energy = (np.sqrt(mom**2 + (vc.mass * vc.amu)**2) - (vc.mass * vc.amu)) / vc.mass
# beta based on brho and q/m
beta = ((vc.c * vc.z * brho) / (vc.amu * vc.mass)) / np.sqrt(((vc.c * vc.z * brho) /
                                                             (vc.amu * vc.mass)) * ((vc.c * vc.z * brho) / (vc.amu * vc.mass)) + 1)
# Lorentz factor
gamma = 1 / np.sqrt(1 - beta**2)
# velocity of the beam [mm/ns]
v = vc.c * beta * 1E-6

print('*******    Initial conditions    *******')
print(f'Proton number of the beam is {vc.z}')
print(f'Mass of the beam is {vc.mass}')
print(f'Mass to charge ratio of the beam is {MoQ}')
print(f'Beam energy is {energy * 1E6} MeV/u')

# ===== Runge-Kutta part =====
# initialize variables
vx = beta * np.sin(a_init / 1000)
vy = beta * np.cos(a_init / 1000)
x = vc.x0 / 1000   # convert mm to m
y = vc.y0 / 1000   # convert mm to m
t = 0.0

# for plot
data_points = []

# for while loop. path is a flight length in mm
stop_flag = 0

print('*******   Runge-Kutta parameters   *******')
print(f'Initial x is {x * 1000} mm')
print(f'Initial y is {y * 1000} mm')

# Precompute constants
step_factor_1 = vc.step_time * vc.c * 1E-9  # around y=0, step time is 100 fs
step_factor_2 = vc.step_time * 100 * vc.c * 1E-9  # other y, step time is 10 ps

while stop_flag < 1:
    # calculate step time of Runge-Kutta
    step = step_factor_1 if ((-0.005 < y) or (y < 0.005)) and (x > 0) else step_factor_2

    # calculate k1
    x_k1 = step * Diff.diff_x(vx)
    y_k1 = step * Diff.diff_y(vy)
    vx_k1 = step * Diff.diff_vx(x + x_k1, y + y_k1, vy, gamma)
    vy_k1 = step * Diff.diff_vy(x + x_k1, y + y_k1, vx, gamma)

    # calculate k2
    x_k2 = step * Diff.diff_x(vx + vx_k1 / 2)
    y_k2 = step * Diff.diff_y(vy + vy_k1 / 2)
    vx_k2 = step * Diff.diff_vx(x + x_k2 / 2, y + y_k2 / 2, vy + vy_k1 / 2, gamma)
    vy_k2 = step * Diff.diff_vy(x + x_k2 / 2, y + y_k2 / 2, vx + vx_k1 / 2, gamma)

    # calculate k3
    x_k3 = step * Diff.diff_x(vx + vx_k2 / 2)
    y_k3 = step * Diff.diff_y(vy + vy_k2 / 2)
    vx_k3 = step * Diff.diff_vx(x + x_k3 / 2, y + y_k3 / 2, vy + vy_k2 / 2, gamma)
    vy_k3 = step * Diff.diff_vy(x + x_k3 / 2, y + y_k3 / 2, vx + vx_k2 / 2, gamma)

    # calculate k4
    x_k4 = step * Diff.diff_x(vx + vx_k3)
    y_k4 = step * Diff.diff_y(vy + vy_k3)
    vx_k4 = step * Diff.diff_vx(x + x_k4, y + y_k4, vy + vy_k3, gamma)
    vy_k4 = step * Diff.diff_vy(x + x_k4, y + y_k4, vx + vx_k3, gamma)

    # calculate next steps
    vx_next = vx + (vx_k1 + 2 * vx_k2 + 2 * vx_k3 + vx_k4) / 6
    vy_next = vy + (vy_k1 + 2 * vy_k2 + 2 * vy_k3 + vy_k4) / 6
    x_next = x + (x_k1 + 2 * x_k2 + 2 * x_k3 + x_k4) / 6
    y_next = y + (y_k1 + 2 * y_k2 + 2 * y_k3 + y_k4) / 6

    data_points.append([t*1E9/vc.c, x * 1000, y * 1000])

    # update variables
    t += step

    if (y < 0) and (x > 0):
        stop_flag = 0.5

    vx = vx_next
    vy = vy_next
    x = x_next
    y = y_next

    if (stop_flag == 0.5) and (y > 0):
        stop_flag = 1

print('*******   Revolution time   *******')
print(f'{t*1E9/vc.c} ns')

# # Output to file
# np.savetxt('main_output.dat', data_points, fmt='%.6f', header='Time(ms) X(mm) Y(mm)', comments='')

# Plot data
plot_data = np.array(data_points)
plot_t, plot_x, plot_y = plot_data[:, 0], plot_data[:, 1], plot_data[:, 2]

print('*******   The number of iterations   *******')
print(f'{len(plot_t)} times')

box2 = plt.figure(figsize=(17.5, 5))
fig1 = box2.add_subplot(1, 3, 1)
plt.plot(plot_x, plot_y, linewidth=1, label='x vs y')
plt.legend()
fig2 = box2.add_subplot(1, 3, 2)
plt.plot(plot_t, plot_x, linewidth=1, label='t vs x')
plt.legend()
fig3 = box2.add_subplot(1, 3, 3)
plt.plot(plot_t, plot_y, linewidth=1, label='t vs y')
plt.legend()

# # comparing yano-results and abe-san-results
# box3 = plt.figure(figsize=(7, 6))
# fig3_1 = box3.add_subplot(1, 1, 1)
# fig3_1 = plt.plot(plot_x, plot_y, linewidth=1, label='yano')
# plt.legend()
# # open and read abe-san-result
# openfile = open('./search_output/kidou_emi_0_dp_0.00_mx_0.dat', 'rt')
# x_list = []
# y_list = []
# for line in openfile:
#     data = line[:-1].split(' ')
#     x_list.append(float(data[1]))
#     y_list.append(float(data[2]))
# fig3_1 = plt.plot(x_list, y_list, color='RED', linewidth=1, label='abe-san')
# magnet_pos_x_flat = []
# for row in vp.magnet_pos_x:
#     magnet_pos_x_flat.extend(row)
# magnet_pos_y_flat = []
# for row in vp.magnet_pos_y:
#     magnet_pos_y_flat.extend(row)
# fig3_1 = plt.plot(magnet_pos_x_flat, magnet_pos_y_flat, 'x',
#                   color='GREEN', label='magnets')
# fig3_1 = plt.plot()
# plt.legend()

# Execution time measurement is stopped
end = time.perf_counter()
print('*******   Execution time is   *******')
print((end - start))

plt.show()
