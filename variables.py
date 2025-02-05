import numpy as np

# *******   magnets   *******
# current of the main magnets
# main_current = 1918.1205
# main_current = 1918.4094  # this value is the best fit with abe-san's code
# main_current = 1911.5 # most closed orbit of z-motion
main_current = 1910.1 # most closed orbit of x-motion

# calculation is normalized by this current
main_base = 1915.0

# set magnetic flux density
# B0z = 1.182166
current_ratio = main_current/main_base

# the positions of dipole magnets
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

magnet_pos_z = ((0, 0, 0, 0),
                (0, 0, 0, 0),
                (0, 0, 0, 0),
                (0, 0, 0, 0),
                (0, 0, 0, 0),
                (0, 0, 0, 0))

angles = np.arange(6) * 60
rad1 = 7.5 + angles
rad2 = 22.5 + angles
rad3 = 37.5 + angles
rad4 = 52.5 + angles

bend_angle = tuple(zip(rad1, rad2, rad3, rad4))