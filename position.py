# import some libraries
import numpy as np

# import variables defined by myself
import variables_position as vp
import variables_conditions as vc


def Position(x, y):
    # for sector1
    x_diff1 = np.zeros(4)
    y_diff1 = np.zeros(4)
    x_particle1 = np.zeros(4)
    y_particle1 = np.zeros(4)
    for i in range(4):
        # x_diff is the difference between input-x and the magnet's central x
        x_diff1[i] = x-vp.magnet_pos_x[0][i]
        # y_diff is the difference between input-y and the magnet's central y
        y_diff1[i] = y-vp.magnet_pos_y[0][i]
        # koko kaeru
        x_particle1[i] = x_diff1[i]*np.cos(vp.bending_angle[0][i]) - \
            y_diff1[i]*np.sin(-1*vp.bending_angle[0][i])
        y_particle1[i] = x_diff1[i]*np.sin(vp.bending_angle[0][i]) + \
            y_diff1[i]*np.cos(-1*vp.bending_angle[0][i])

    # for sector2
    x_diff2 = np.zeros(4)
    y_diff2 = np.zeros(4)
    for i in range(4):
        x_diff2[i] = x-vp.magnet_pos_x[1][i]
        y_diff2[i] = y-vp.magnet_pos_y[1][i]

    x_particle2 = np.zeros(4)
    y_particle2 = np.zeros(4)
    for i in range(4):
        x_particle2[i] = x_diff2[i]*np.cos(vp.bending_angle[1][i]) - \
            y_diff2[i]*np.sin(-1*vp.bending_angle[1][i])
        y_particle2[i] = x_diff2[i]*np.sin(vp.bending_angle[1][i]) + \
            y_diff2[i]*np.cos(-1*vp.bending_angle[1][i])

    # for sector3
    x_diff3 = np.zeros(4)
    y_diff3 = np.zeros(4)
    for i in range(4):
        x_diff3[i] = x-vp.magnet_pos_x[2][i]
        y_diff3[i] = y-vp.magnet_pos_y[2][i]

    x_particle3 = np.zeros(4)
    y_particle3 = np.zeros(4)
    for i in range(4):
        x_particle3[i] = x_diff3[i]*np.cos(vp.bending_angle[2][i]) - \
            y_diff3[i]*np.sin(-1*vp.bending_angle[2][i])
        y_particle3[i] = x_diff3[i]*np.sin(vp.bending_angle[2][i]) + \
            y_diff3[i]*np.cos(-1*vp.bending_angle[2][i])

    # for sector4
    x_diff4 = np.zeros(4)
    y_diff4 = np.zeros(4)
    for i in range(4):
        x_diff4[i] = x-vp.magnet_pos_x[3][i]
        y_diff4[i] = y-vp.magnet_pos_y[3][i]

    x_particle4 = np.zeros(4)
    y_particle4 = np.zeros(4)
    for i in range(4):
        x_particle4[i] = x_diff4[i]*np.cos(vp.bending_angle[3][i]) - \
            y_diff4[i]*np.sin(-1*vp.bending_angle[3][i])
        y_particle4[i] = x_diff4[i]*np.sin(vp.bending_angle[3][i]) + \
            y_diff4[i]*np.cos(-1*vp.bending_angle[3][i])

    # for sector5
    x_diff5 = np.zeros(4)
    y_diff5 = np.zeros(4)
    for i in range(4):
        x_diff5[i] = x-vp.magnet_pos_x[4][i]
        y_diff5[i] = y-vp.magnet_pos_y[4][i]

    x_particle5 = np.zeros(4)
    y_particle5 = np.zeros(4)
    for i in range(4):
        x_particle5[i] = x_diff5[i]*np.cos(vp.bending_angle[4][i]) - \
            y_diff5[i]*np.sin(-1*vp.bending_angle[4][i])
        y_particle5[i] = x_diff5[i]*np.sin(vp.bending_angle[4][i]) + \
            y_diff5[i]*np.cos(-1*vp.bending_angle[4][i])

    # for sector6
    x_diff6 = np.zeros(4)
    y_diff6 = np.zeros(4)
    for i in range(4):
        x_diff6[i] = x-vp.magnet_pos_x[5][i]
        y_diff6[i] = y-vp.magnet_pos_y[5][i]

    x_particle6 = np.zeros(4)
    y_particle6 = np.zeros(4)
    for i in range(4):
        x_particle6[i] = x_diff6[i]*np.cos(vp.bending_angle[5][i]) - \
            y_diff6[i]*np.sin(-1*vp.bending_angle[5][i])
        y_particle6[i] = x_diff6[i]*np.sin(vp.bending_angle[5][i]) + \
            y_diff6[i]*np.cos(-1*vp.bending_angle[5][i])

    pos = np.zeros((6, 2, 4))
    for i in range(4):
        pos[0][0][i] = x_particle1[i]
        pos[0][1][i] = y_particle1[i]

        pos[1][0][i] = x_particle2[i]
        pos[1][1][i] = y_particle2[i]

        pos[2][0][i] = x_particle3[i]
        pos[2][1][i] = y_particle3[i]

        pos[3][0][i] = x_particle4[i]
        pos[3][1][i] = y_particle4[i]

        pos[4][0][i] = x_particle5[i]
        pos[4][1][i] = y_particle5[i]

        pos[5][0][i] = x_particle6[i]
        pos[5][1][i] = y_particle6[i]

    return pos


pos = Position(vc.x0, vc.y0)
print(pos)
