import numpy as np

# position of the center of the magnets
# 4 magnets x 6 sectors

# magnet_pos_x is the center of magnet's coordinate in the X-axis
magnet_pos_x = []
# magnet_pos_y is the center of magnet's coordinate in the Y-axis
magnet_pos_y = []

magnet_pos_x = ((9232.420563, 8805.253915, 7980.153668, 6812.990034),
                (2419.652338, 823.734813, -824.801225, -2419.208113),
                (-6812.546113, -7979.587097, -8804.80969, -9231.976338),
                (-9231.976338, -8804.809690, -7979.587097, -6812.546113),
                (-2419.208113, -825.000481, 825.444706, 2419.652338),
                (6812.990338, 7980.031322, 8805.253915, 9232.420563))

magnet_pos_y = ((2536.110167, 4130.317798, 5559.485814, 6726.686242),
                (9263.181119, 9677.359839, 9690.373999, 9263.181119),
                (6726.686242, 5559.645258, 4130.317798, 2536.10167),
                (-2536.879587, -4131.087218, -5560.414678, -6727.455662),
                (-9263.950539, -9691.117187, -9691.117187, -9263.950539),
                (-6727.455662, -5560.414678, -4131.087218, -2536.879587))

# bending angle of the magnets in the X-Y plane
rad1 = np.zeros(6)
rad2 = np.zeros(6)
rad3 = np.zeros(6)
rad4 = np.zeros(6)

for i in range(6):
    rad1[i] = (7.5+i*60)
    rad2[i] = (22.5+i*60)
    rad3[i] = (37.5+i*60)
    rad4[i] = (52.5+i*60)


bend_angle = []

bend_angle = ((rad1[0], rad2[0], rad3[0], rad4[0]),
              (rad1[1], rad2[1], rad3[1], rad4[1]),
              (rad1[2], rad2[2], rad3[2], rad4[2]),
              (rad1[3], rad2[3], rad3[3], rad4[3]),
              (rad1[4], rad2[4], rad3[4], rad4[4]),
              (rad1[5], rad2[5], rad3[5], rad4[5]))

