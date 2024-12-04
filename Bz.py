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

# x_range=np.arange(-10000,10000,100)
# y_range=np.arange(-10000,10000,100)
# x,y=np.meshgrid(x_range,y_range)
# BforXplane_vectorized = np.vectorize(BforXplane)  # ベクトル化
# z = BforXplane_vectorized(x, y)

# plt.figure(figsize=(10.5, 9))
# plt.contourf(x, y, z, levels=100, cmap='viridis')  # 等高線プロット
# plt.colorbar(label='B(x, y)')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Mapping of B(x, y)')
# plt.show()
  
# fig = plt.figure()
# ax=fig.add_subplot(projection='3d')
# ax.plot_wireframe(x,y,z)
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("B")
# plt.show()