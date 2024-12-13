# import some libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# import variables defined by myself
import variables_conditions as vc
import variables_Bz as vb
import variables_position as vp

# import functions made by myself
import FourGaussian
import Enge
import position


# if magnet has trim coils, Bz : B0z * B-X curve * trim * fringe
# else if there is no trim coils, Bz : B0z * B-X curve * fringe
# else input parameter is wrong...
def bz(trim, x_diff, x, y):
    # print(x,y)
    # if (abs(x) < 0.2*1e3) and (abs(y) < 1000*1e3):
    if (abs(x) < 0.2) and (abs(y) < 1000):
        if trim == 'trim':
            Bz = vc.B0z *\
                (vb.para_twz[0]+vb.para_twz[1]*(x) +
                 vb.para_twz[2]*(x)**2+vb.para_twz[3]*(x)**3 +
                 vb.para_twz[4]*(x)**4+vb.para_twz[5]*(x)**5 +
                 vb.para_twz[6]*(x)**6+vb.para_twz[7]*(x)**7 +
                 vb.para_twz[8]*(x)**8)
            # if Bz < 0:
            #     print('Im here')
            #     return 0
            # else:
            Bz = Bz+FourGaussian.Btrim_Sum(x_diff)
            Bz = Bz*Enge.enge(y)*-1
            return Bz
        elif trim == 'no_trim':
            Bz = vc.B0z *\
                (vb.para_twz[0]+vb.para_twz[1]*(x) +
                 vb.para_twz[2]*(x)**2+vb.para_twz[3]*(x)**3 +
                 vb.para_twz[4]*(x)**4+vb.para_twz[5]*(x)**5 +
                 vb.para_twz[6]*(x)**6+vb.para_twz[7]*(x)**7 +
                 vb.para_twz[8]*(x)**8)
            # if Bz < 0:
            #     return 0
            # else:
            Bz = Bz*Enge.enge(y)*-1
            return Bz
        else:
            print('You put a wrong input parameter!!')
            print('Check Bz.py')
            return 0
    else:
        return 0


def BforXplane(x, y):
    # Prepare required variables
    magnet_x = np.array(vp.magnet_pos_x)
    magnet_y = np.array(vp.magnet_pos_y)
    bending_angle = np.array(vp.bend_angle)
    
    sectors = []
    BforXplane = 0.0
    
    # Convert x and y to mm from m
    x_mm, y_mm = x * 1000, y * 1000

    # Calculate positions for each magnet sector
    for i in range(magnet_x.shape[0]):  # Loop over sectors
        sector = []
        for j in range(magnet_x.shape[1]):  # Loop over magnets in a sector
            sector.append(
                position.Position(
                    x_mm, y_mm, magnet_x[i, j], magnet_y[i, j], bending_angle[i, j]
                )
            )
        sectors.append(sector)

    # Calculate Bz for each sector and magnet
    for i, sector in enumerate(sectors):
        for j, magnet in enumerate(sector):
            trim = 'trim' if j % 4 == 0 or j % 4 == 3 else 'no_trim'
            BforXplane += bz(trim, magnet[0][0], magnet[0][1], magnet[0][2])

    return BforXplane



# # for plotting the map of magnetic field
# x_range = np.arange(4.5, 11.5, 0.25)  # unit is m
# y_range = np.arange(1, 8, 0.25)  # unit is m
# X, Y = np.meshgrid(x_range, y_range)
# Z = np.zeros_like(X)
# # calculate BforXplane for each x_range and y_range
# for i in range(X.shape[0]):
#     for j in range(X.shape[1]):
#         x, y = X[i, j], Y[i, j]
#         Z[i, j] = BforXplane(x, y)
# # 3D plot
# fig = plt.figure(figsize=(10, 7))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_wireframe(X*1000, Y*1000, -Z, color='blue', linewidth=0.5, rcount=200, ccount=200)
# # label and title
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('BforXplane')
# plt.show()


# # for plotting bz()
# plot_x = np.arange(-0.25, 0.25, 0.001)
# plot_y = np.arange(-250, 250, 1)

# # put 0 into plot_B if B is smaller than 0
# condition = vc.B0z * (vb.para_twz[0] + vb.para_twz[1] * (plot_x-0) +
#                       vb.para_twz[2] * ((plot_x-0)**2) + vb.para_twz[3] * ((plot_x-0)**3) +
#                       vb.para_twz[4] * ((plot_x-0)**4) + vb.para_twz[5] * ((plot_x-0)**5) +
#                       vb.para_twz[6] * ((plot_x-0)**6) + vb.para_twz[7] * ((plot_x-0)**7) +
#                       vb.para_twz[8] * ((plot_x-0)**8) +
#                       FourGaussian.Btrim_Sum((plot_x-0)*1000)) * Enge.enge(plot_y)*-1< 0

# plot_B = np.where(
#     condition,
#     0,  # when condition is TRUE
#     vc.B0z * (vb.para_twz[0] + vb.para_twz[1] * (plot_x-0) +
#               vb.para_twz[2] * ((plot_x-0)**2) + vb.para_twz[3] * ((plot_x-0)**3) +
#               vb.para_twz[4] * ((plot_x-0)**4) + vb.para_twz[5] * ((plot_x-0)**5) +
#               vb.para_twz[6] * ((plot_x-0)**6) + vb.para_twz[7] * ((plot_x-0)**7) +
#               vb.para_twz[8] * ((plot_x-0)**8)) +
#               FourGaussian.Btrim_Sum((plot_x-0)*1000) * Enge.enge(plot_y)*-1   # when condition is FALSE
# )

# fig2 = plt.figure(figsize=(10.5, 9))
# bx = fig2.add_subplot(1, 1, 1)
# bx.plot(plot_x*1000, plot_B, label='x-axis')
# bx.plot(plot_y, plot_B, label='y-axis')
# bx.legend()
# plt.show()
