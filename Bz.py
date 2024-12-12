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
    if (abs(x) < 0.2) and (abs(y) < 1000):
        if trim == 'trim':
            Bz = vc.B0z *\
                (vb.para_twz[0]+vb.para_twz[1]*(x) +
                 vb.para_twz[2]*(x)**2+vb.para_twz[3]*(x)**3 +
                 vb.para_twz[4]*(x)**4+vb.para_twz[5]*(x)**5 +
                 vb.para_twz[6]*(x)**6+vb.para_twz[7]*(x)**7 +
                 vb.para_twz[8]*(x)**8)
            if Bz < 0:
                return 0
            else:
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
            if Bz < 0:
                return 0
            else:
                Bz = Bz*Enge.enge(y)*-1
                return Bz
        else:
            print('You put a wrong input parameter!!')
            print('Check Bz.py')
            return 0
    else:
        return 0


def BforXplane(x, y):
    magnet_x = vp.magnet_pos_x
    magnet_y = vp.magnet_pos_y
    bending_angle = vp.bend_angle
    # print(magnet_x[0][0])
    # print(magnet_y[0][0])
    # print(bending_angle[0][0])

    # calculate coordinate for each bending magnet
    # sector 1
    # note that x and y are converted to mm from m
    sector1_1 = position.Position(
        x*1000, y*1000, magnet_x[0][0], magnet_y[0][0], bending_angle[0][0])
    sector1_2 = position.Position(
        x*1000, y*1000, magnet_x[0][1], magnet_y[0][1], bending_angle[0][1])
    sector1_3 = position.Position(
        x*1000, y*1000, magnet_x[0][2], magnet_y[0][2], bending_angle[0][2])
    sector1_4 = position.Position(
        x*1000, y*1000, magnet_x[0][3], magnet_y[0][3], bending_angle[0][3])
    # sector 2
    sector2_1 = position.Position(
        x*1000, y*1000, magnet_x[1][0], magnet_y[1][0], bending_angle[1][0])
    sector2_2 = position.Position(
        x*1000, y*1000, magnet_x[1][1], magnet_y[1][1], bending_angle[1][1])
    sector2_3 = position.Position(
        x*1000, y*1000, magnet_x[1][2], magnet_y[1][2], bending_angle[1][2])
    sector2_4 = position.Position(
        x*1000, y*1000, magnet_x[1][3], magnet_y[1][3], bending_angle[1][3])
    # sector 3
    sector3_1 = position.Position(
        x*1000, y*1000, magnet_x[2][0], magnet_y[2][0], bending_angle[2][0])
    sector3_2 = position.Position(
        x*1000, y*1000, magnet_x[2][1], magnet_y[2][1], bending_angle[2][1])
    sector3_3 = position.Position(
        x*1000, y*1000, magnet_x[2][2], magnet_y[2][2], bending_angle[2][2])
    sector3_4 = position.Position(
        x*1000, y*1000, magnet_x[2][3], magnet_y[2][3], bending_angle[2][3])
    # sector 4
    sector4_1 = position.Position(
        x*1000, y*1000, magnet_x[3][0], magnet_y[3][0], bending_angle[3][0])
    sector4_2 = position.Position(
        x*1000, y*1000, magnet_x[3][1], magnet_y[3][1], bending_angle[3][1])
    sector4_3 = position.Position(
        x*1000, y*1000, magnet_x[3][2], magnet_y[3][2], bending_angle[3][2])
    sector4_4 = position.Position(
        x*1000, y*1000, magnet_x[3][3], magnet_y[3][3], bending_angle[3][3])
    # sector 5
    sector5_1 = position.Position(
        x*1000, y*1000, magnet_x[4][0], magnet_y[4][0], bending_angle[4][0])
    sector5_2 = position.Position(
        x*1000, y*1000, magnet_x[4][1], magnet_y[4][1], bending_angle[4][1])
    sector5_3 = position.Position(
        x*1000, y*1000, magnet_x[4][2], magnet_y[4][2], bending_angle[4][2])
    sector5_4 = position.Position(
        x*1000, y*1000, magnet_x[4][3], magnet_y[4][3], bending_angle[4][3])
    # sector 6
    sector6_1 = position.Position(
        x*1000, y*1000, magnet_x[5][0], magnet_y[5][0], bending_angle[5][0])
    sector6_2 = position.Position(
        x*1000, y*1000, magnet_x[5][1], magnet_y[5][1], bending_angle[5][1])
    sector6_3 = position.Position(
        x*1000, y*1000, magnet_x[5][2], magnet_y[5][2], bending_angle[5][2])
    sector6_4 = position.Position(
        x*1000, y*1000, magnet_x[5][3], magnet_y[5][3], bending_angle[5][3])

    # [0][0] means x2 in abe-san's code, [0][1] means xpos in abe-san's code,
    # and [0][2] means dy in abe-san's code
    BforXplane = bz('trim', sector1_1[0][0], sector1_1[0][1], sector1_1[0][2]) + \
        bz('no_trim', sector1_2[0][0], sector1_2[0][1], sector1_2[0][2]) + \
        bz('no_trim', sector1_3[0][0], sector1_3[0][1], sector1_3[0][2]) + \
        bz('trim', sector1_4[0][0], sector1_4[0][1], sector1_4[0][2])
    BforXplane = BforXplane + \
        bz('trim', sector2_1[0][0], sector2_1[0][1], sector2_1[0][2]) + \
        bz('no_trim', sector2_2[0][0], sector2_2[0][1], sector2_2[0][2]) + \
        bz('no_trim', sector2_3[0][0], sector2_3[0][1], sector2_3[0][2]) + \
        bz('trim', sector2_4[0][0], sector2_4[0][1], sector2_4[0][2])
    BforXplane = BforXplane + \
        bz('trim', sector3_1[0][0], sector3_1[0][1], sector3_1[0][2]) + \
        bz('no_trim', sector3_2[0][0], sector3_2[0][1], sector3_2[0][2]) + \
        bz('no_trim', sector3_3[0][0], sector3_3[0][1], sector3_3[0][2]) + \
        bz('trim', sector3_4[0][0], sector3_4[0][1], sector3_4[0][2])
    BforXplane = BforXplane + \
        bz('trim', sector4_1[0][0], sector4_1[0][1], sector4_1[0][2]) + \
        bz('no_trim', sector4_2[0][0], sector4_2[0][1], sector4_2[0][2]) + \
        bz('no_trim', sector4_3[0][0], sector4_3[0][1], sector4_3[0][2]) + \
        bz('trim', sector4_4[0][0], sector4_4[0][1], sector4_4[0][2])
    BforXplane = BforXplane + \
        bz('trim', sector5_1[0][0], sector5_1[0][1], sector5_1[0][2]) + \
        bz('no_trim', sector5_2[0][0], sector5_2[0][1], sector5_2[0][2]) + \
        bz('no_trim', sector5_3[0][0], sector5_3[0][1], sector5_3[0][2]) + \
        bz('trim', sector5_4[0][0], sector5_4[0][1], sector5_4[0][2])
    BforXplane = BforXplane + \
        bz('trim', sector6_1[0][0], sector6_1[0][1], sector6_1[0][2]) + \
        bz('no_trim', sector6_2[0][0], sector6_2[0][1], sector6_2[0][2]) + \
        bz('no_trim', sector6_3[0][0], sector6_3[0][1], sector6_3[0][2]) + \
        bz('trim', sector6_4[0][0], sector6_4[0][1], sector6_4[0][2])

    return BforXplane


# # for plotting the map of magnetic field
# x_range = np.arange(4.5, 11.5, 0.005)  # unit is m
# y_range = np.arange(1, 8, 0.005)  # unit is m
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
# ax.plot_wireframe(X*1000, Y*1000, -Z, color='blue', linewidth=0.5, rcount=100, ccount=100)
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
