import numpy as np
import matplotlib.pyplot as plt

import variables_conditions as vc

import Diff

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

# m/q
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
px = mom
py = 0.0
x = vc.x0
y = vc.y0
t = 0.0
# for plot
plot_x = []
plot_y = []
plot_t = []
# for wile loop
path = 0.0
number_step = 0.0

while path < 6000:
    # calculate k1
    px_k1 = vc.step_time*Diff.diff_px(x, y, v, gamma, 'trim')
    py_k1 = vc.step_time*Diff.diff_px(x, y, v, gamma, 'trim')
    x_k1 = vc.step_time*Diff.diff_x(px)
    y_k1 = vc.step_time*Diff.diff_x(py)
    # calculate k2
    px_k2 = vc.step_time*Diff.diff_px(x+x_k1/2, y+y_k1/2, v, gamma, 'trim')
    py_k2 = vc.step_time*Diff.diff_px(x+x_k1/2, y+y_k1/2, v, gamma, 'trim')
    x_k2 = vc.step_time*Diff.diff_x(px+px_k1/2)
    y_k2 = vc.step_time*Diff.diff_x(py+py_k1/2)
    # calculate k3
    px_k3 = vc.step_time*Diff.diff_px(x+x_k2/2, y+y_k2/2, v, gamma, 'trim')
    py_k3 = vc.step_time*Diff.diff_px(x+x_k2/2, y+y_k2/2, v, gamma, 'trim')
    x_k3 = vc.step_time*Diff.diff_x(px+px_k2/2)
    y_k3 = vc.step_time*Diff.diff_x(py+py_k2/2)
    # calculate k4
    px_k4 = vc.step_time*Diff.diff_px(x+x_k3, y+y_k3, v, gamma, 'trim')
    py_k4 = vc.step_time*Diff.diff_px(x+x_k3, y+y_k3, v, gamma, 'trim')
    x_k4 = vc.step_time*Diff.diff_x(px+px_k3)
    y_k4 = vc.step_time*Diff.diff_x(py+py_k3)
    # calculate next steps
    px_next = px+(px_k1 + 2*px_k2 + 2*px_k3 + px_k4)/6
    py_next = py+(py_k1 + 2*py_k2 + 2*py_k3 + py_k4)/6
    x_next = x+(x_k1 + 2*x_k2 + 2*x_k3 + x_k4)/6
    y_next = y+(y_k1 + 2*y_k2 + 2*y_k3 + y_k4)/6
    plot_x.append(x)
    plot_y.append(y)
    plot_t.append(t)
    # update variables
    px = px_next
    py = py_next
    x = x_next
    y = y_next
    t += vc.step_time
    number_step += 1
    path += number_step*v

box = plt.figure(figsize=(14, 6))
fig1 = box.add_subplot(1, 2, 1)
fig1 = plt.plot(plot_x, plot_y, linewidth=1, label='qx vs qy')
plt.legend()
fig2 = box.add_subplot(1, 2, 2)
fig2 = plt.plot(plot_t, plot_x, linewidth=1, label='t vs qx')
plt.legend()
plt.show()
