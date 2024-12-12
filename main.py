import time
# for measuring execution time
start = time.perf_counter()

import position
import Bz
import Diff
import variables_conditions as vc
import matplotlib.pyplot as plt
import variables_position as vp
import numpy as np


#  set brho
j = 0.0
dp = 0.1*j
brho = vc.brho0*(1+dp/100)

# set emittance
l = 0.0
emi = l*20   # 20 steps
a_max = np.sqrt(emi/vc.betax)*0.95
a_bin = 2*a_max/15
k = 0.0
a_init = k*a_bin - a_max

# set dispersion [m/%] in the ring
dispersion = -3.3494*dp**6 + 2.2088*dp**5 + 1.5699*dp**4 - \
    0.731*dp**3 - 0.1573*dp**2 + 0.0121*dp + 7.0737

# q/m
MoQ = vc.mass/vc.z
# momentum based on brho
mom = vc.c*vc.z*brho*1E-6
# beam energy based on momentum
energy = (np.sqrt(mom**2+(vc.mass*vc.amu)**2)-(vc.mass*vc.amu))/vc.mass
# beta based on brho and q/m
beta = ((vc.c*vc.z*brho)/(vc.amu*vc.mass))/np.sqrt(((vc.c*vc.z*brho) /
                                                    (vc.amu*vc.mass))*((vc.c*vc.z*brho)/(vc.amu*vc.mass))+1)
# Lorentz factor
gamma = 1/np.sqrt(1-beta**2)
# velocity of the beam [mm/ns]
v = vc.c*beta*1E-6

print('*******    Initial conditions    *******')
print(f'Proton number of the beam is {vc.z}')
print(f'Mass of the beam is {vc.mass}')
print(f'Mass to charge ratio of the beam is {MoQ}')
print(f'Beam energy is {energy*1E6} MeV/u')

# ===== Runge-Kutta part =====
# initialize variables
vx = beta*np.sin(a_init/1000)
vy = beta*np.cos(a_init/1000)
x = vc.x0/1000   # convert mm to m
y = vc.y0/1000   # convert mm to m
t = 0.0
step = vc.step_time*100*vc.c*1E-9
# for plot
plot_x = []
plot_y = []
plot_vx = []
plot_vy = []
plot_t = []
# for wile loop. path is a flight length in mm
path = 0.0

print('*******   Runge-Kutta parameters   *******')
print(f'Initial vx is {vx}')
print(f'Initial vy is {vy}')
print(f'Initial x is {x}')
print(f'Initial y is {y}')

while path < 603.5:
    # calculate k1
    x_k1 = step*Diff.diff_x(vx)
    y_k1 = step*Diff.diff_y(vy)
    vx_k1 = step*Diff.diff_vx(x+x_k1, y+y_k1, vy, gamma)
    vy_k1 = step*Diff.diff_vy(x+x_k1, y+y_k1, vx, gamma)
    # calculate k2
    x_k2 = step*Diff.diff_x(vx+vx_k1/2)
    y_k2 = step*Diff.diff_y(vy+vy_k1/2)
    vx_k2 = step*Diff.diff_vx(x+x_k2/2, y+y_k2/2, vy+vy_k1/2, gamma)
    vy_k2 = step*Diff.diff_vy(x+x_k2/2, y+y_k2/2, vx+vx_k1/2, gamma)
    # calculate k3
    x_k3 = step*Diff.diff_x(vx+vx_k2/2)
    y_k3 = step*Diff.diff_y(vy+vy_k2/2)
    vx_k3 = step*Diff.diff_vx(x+x_k3/2, y+y_k3/2, vy+vy_k2/2, gamma)
    vy_k3 = step*Diff.diff_vy(x+x_k3/2, y+y_k3/2, vx+vx_k2/2, gamma)
    # calculate k4
    x_k4 = step*Diff.diff_x(vx+vx_k3)
    y_k4 = step*Diff.diff_y(vy+vy_k3)
    vx_k4 = step*Diff.diff_vx(x+x_k4, y+y_k4, vy+vy_k3, gamma)
    vy_k4 = step*Diff.diff_vy(x+x_k4, y+y_k4, vx+vx_k3, gamma)
    # calculate next steps
    vx_next = vx+(vx_k1 + 2*vx_k2 + 2*vx_k3 + vx_k4)/6
    vy_next = vy+(vy_k1 + 2*vy_k2 + 2*vy_k3 + vy_k4)/6
    x_next = x+(x_k1 + 2*x_k2 + 2*x_k3 + x_k4)/6
    y_next = y+(y_k1 + 2*y_k2 + 2*y_k3 + y_k4)/6
    plot_x.append(x*1000)
    plot_y.append(y*1000)
    plot_vx.append(vx)
    plot_vy.append(vy)
    plot_t.append(t*1000)
    # update variables
    path += v*vc.step_time
    t += vc.step_time
    vx = vx_next
    vy = vy_next
    x = x_next
    y = y_next

print('*******   Revolution time   *******')
print(f'{t*1000} ns')

# making output file
outfile = open('main_outpit.dat', 'w')
for i in range(len(plot_x)):
    outfile.write(f'{plot_t[i]} {plot_x[i]} {plot_y[i]}')
outfile.close

box2 = plt.figure(figsize=(17.5, 5))
fig1 = box2.add_subplot(1, 3, 1)
fig1 = plt.plot(plot_x, plot_y, linewidth=1, label='x vs y')
plt.legend()
fig2 = box2.add_subplot(1, 3, 2)
fig2 = plt.plot(plot_t, plot_x, linewidth=1, label='t vs x')
plt.legend()
fig3 = box2.add_subplot(1, 3, 3)
fig2 = plt.plot(plot_t, plot_y, linewidth=1, label='t vs y')
plt.legend()

# comparing yano-results and abe-san-results
box3 = plt.figure(figsize=(7, 6))
fig3_1 = box3.add_subplot(1, 1, 1)
fig3_1 = plt.plot(plot_x, plot_y, linewidth=1, label='yano')
plt.legend()
# open and read abe-san-result
openfile = open('./search_output/kidou_emi_0_dp_0.00_mx_0.dat', 'rt')
x_list = []
y_list = []
for line in openfile:
    data = line[:-1].split(' ')
    x_list.append(float(data[1]))
    y_list.append(float(data[2]))

fig3_1 = plt.plot(x_list, y_list, color='RED', linewidth=1, label='abe-san')
magnet_pos_x_flat = []
for row in vp.magnet_pos_x:
    magnet_pos_x_flat.extend(row)
magnet_pos_y_flat = []
for row in vp.magnet_pos_y:
    magnet_pos_y_flat.extend(row)
fig3_1 = plt.plot(magnet_pos_x_flat, magnet_pos_y_flat, 'x',
                  color='GREEN', label='magnets')
fig3_1 = plt.plot()
plt.legend()

box4 = plt.figure(figsize=(17.5, 5))
fig4_1 = box4.add_subplot(1, 2, 1)
fig4_1 = plt.plot(plot_t, plot_vx, linewidth=1, label='t vs vx')
plt.legend()
fig4_2 = box4.add_subplot(1, 2, 2)
fig4_2 = plt.plot(plot_t, plot_vy, linewidth=1, label='t vs vy')
plt.legend()

# execution time measurement is stopped
end = time.perf_counter()
print('*******   Execution time is   *******')
print((end-start))

plt.show()
