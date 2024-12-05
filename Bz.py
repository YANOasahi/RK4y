# import some libraries
import numpy as np
import matplotlib.pyplot as plt

# import variables defined by myself
import variables_conditions as vc
import variables_Bz as vb

# import functions made by myself
import FourGaussian
import Enge
import position


# if magnet has trim coils, Bz : B0z * B-X curve * trim * fringe
# else if there is no trim coils, Bz : B0z * B-X curve * fringe
# else input parameter is wrong...
def bz(trim, x, y):
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
            Bz = Bz+FourGaussian.Btrim_Sum(x)
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
        Bz = 0
        return Bz


def BforXplane(x, y):
    pos = position.Position(x, y)
    BforXplane = bz('trim', pos[0][0][0], pos[0][1][0]) + \
        bz('no_trim', pos[0][0][1], pos[0][1][1]) + \
        bz('no_trim', pos[0][0][2], pos[0][1][2]) + \
        bz('trim', pos[0][0][3], pos[0][1][3])
    BforXplane = BforXplane + \
        bz('trim', pos[1][0][0], pos[1][1][0]) + \
        bz('no_trim', pos[1][0][1], pos[1][1][1]) + \
        bz('no_trim', pos[1][0][2], pos[1][1][2]) + \
        bz('trim', pos[1][0][3], pos[1][1][3])
    BforXplane = BforXplane + \
        bz('trim', pos[2][0][0], pos[2][1][0]) + \
        bz('no_trim', pos[2][0][1], pos[2][1][1]) + \
        bz('no_trim', pos[2][0][2], pos[2][1][2]) + \
        bz('trim', pos[2][0][3], pos[2][1][3])
    BforXplane = BforXplane + \
        bz('trim', pos[3][0][0], pos[3][1][0]) + \
        bz('no_trim', pos[3][0][1], pos[3][1][1]) + \
        bz('no_trim', pos[3][0][2], pos[3][1][2]) + \
        bz('trim', pos[3][0][3], pos[3][1][3])
    BforXplane = BforXplane + \
        bz('trim', pos[4][0][0], pos[4][1][0]) + \
        bz('no_trim', pos[4][0][1], pos[4][1][1]) + \
        bz('no_trim', pos[4][0][2], pos[4][1][2]) + \
        bz('trim', pos[4][0][3], pos[4][1][3])
    BforXplane = BforXplane + \
        bz('trim', pos[5][0][0], pos[5][1][0]) + \
        bz('no_trim', pos[5][0][1], pos[5][1][1]) + \
        bz('no_trim', pos[5][0][2], pos[5][1][2]) + \
        bz('trim', pos[5][0][3], pos[5][1][3])
    return BforXplane


x_range = np.arange(-10000, 10000, 10)
y_range = np.arange(-10000, 10000, 10)
x, y = np.meshgrid(x_range, y_range)
# Because BforXplane can only accept scalar, we need to vectorize
BforXplane_vectorized = np.vectorize(BforXplane)
z = BforXplane_vectorized(x, y)

fig = plt.figure(figsize=(10.5, 9))
ax = fig.add_subplot(projection='3d')
ax.plot_wireframe(x, y, z)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("B")
plt.show()

# plot_x = np.arange(-0.25, 0.25, 0.001)

# # put 0 into plot_B if B is smaller than 0
# condition = vc.B0z * (vb.para_twz[0] + vb.para_twz[1] * plot_x +
#                       vb.para_twz[2] * (plot_x**2) + vb.para_twz[3] * (plot_x**3) +
#                       vb.para_twz[4] * (plot_x**4) + vb.para_twz[5] * (plot_x**5) +
#                       vb.para_twz[6] * (plot_x**6) + vb.para_twz[7] * (plot_x**7) +
#                       vb.para_twz[8] * (plot_x**8)) < 0

# plot_B = np.where(
#     condition,
#     0,  # when condition is TRUE
#     vc.B0z * (vb.para_twz[0] + vb.para_twz[1] * plot_x +
#               vb.para_twz[2] * (plot_x**2) + vb.para_twz[3] * (plot_x**3) +
#               vb.para_twz[4] * (plot_x**4) + vb.para_twz[5] * (plot_x**5) +
#               vb.para_twz[6] * (plot_x**6) + vb.para_twz[7] * (plot_x**7) +
#               vb.para_twz[8] * (plot_x**8))  # when condition is FALSE
# )

# fig2 = plt.figure(figsize=(10.5, 9))
# bx = fig2.add_subplot(1, 1, 1)
# bx.plot(plot_x, plot_B)
# plt.show()
