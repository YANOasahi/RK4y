import numpy as np

# *******   magnets   *******
# current of the main magnets
# main_current = 1918.1205
main_current = 1918.4094  # this value is the best fit with abe-san's code
# main_current = 1918.383  # return to almost the same position

# current of the trim coils
trim_current = (592.128, 1.958, 91.045, 226.784,
                158.111, 146.1, 117.213, 280.012, 47.809, 279.008)

# calculation is normalized by this current
main_base = 1915.0

# set magnetic flux density
B0z = 1.182166

# parameters for main magnetic field
para_twz = (1*(main_current/main_base), 0.0, 0.0639773,
            0.0, -54.1635, 0.0, 12402.1, 0.0, -651257)

# parameters for fringe field
enge = (0.288396, 1.41928, -0.319549, 0.226168, -0.026047, 0.002317)

# parameters for trim coils
amp = np.array([[-0.00559819, -0.0000670702, -0.000657458, -0.00426336],
                [-0.004793, -0.0000886371, -0.00475243, -0.00252437],
                [-0.002201, -0.001263, -0.0046398, -0.00483284],
                [-0.005188, -0.001419, -0.00244726, -0.00463077],
                [-0.023787, -0.00235, 0.0180907, -0.00533062],
                [0.0019743, 0.0019491, 0.0256026, -0.0193996],
                [0.0097702, -0.000999232, 0.00363535, -0.00586014],
                [0.0037986, 0.0000549654, 0.00239795, 0.0051618],
                [0.0051218, 0.000169799, 0.000811862, 0.00418408],
                [0.0067407, -0.000154026, 0.000703519, -0.000792968]])

mean = np.array([[-126.476, -11.2285, -95.6711, -66.6196],
                 [-137.427, -3.66384, -82.0093, -44.9389],
                 [-157.189, -21.5678, -117.715, -58.3225],
                 [-111.732, -13.3005, -158.686, -45.0715],
                 [-132.117, -6.34961, -129.043, -46.9839],
                 [188.603, 29.0324, 111.263, 112.976],
                 [161.154, 19.2292, 66.6887, 244.062],
                 [168.821, 17.3227, 74.4467, 114.386],
                 [155.302, 34.453, 81.6906, 103.053],
                 [111.382, 614.047, 86.7074, 90.6008]])

sigma = np.array([[30.9606, -13.9661, -14.8655, 23.9585],
                  [29.6729, -13.7996, -29.2736, 21.5422],
                  [25.2112, -18.9073, -32.4186, 32.4758],
                  [36.8752, -19.4676, -26.2377, 32.927],
                  [33.6369, -22.2167, -30.9078, 36.0638],
                  [24.8122, -21.0937, -39.5128, 34.7987],
                  [62.3175, 20.2542, -34.7429, 103.475],
                  [27.9871, -11.5058, -22.043, 32.8979],
                  [30.3747, -17.4388, -16.447, 25.963],
                  [34.6009, -521.179, -15.345, 15.0778]])

offset = np.array([[10.0687, 9.48442, 10.305, 10.1855],
                   [10.1344, 4.08676, 10.3313, 9.67604],
                   [10.2894, 16.9878, 10.7113, 9.93579],
                   [9.96632, 11.0316, 11.1389, 9.66259],
                   [10.0629, 5.98738, 12.9681, 9.5711],
                   [9.23832, -4.80779, 7.609, 9.26375],
                   [9.8653, 3.40861, 10.6119, 6.11265],
                   [9.80403, 5.13617, 10.4978, 9.07634],
                   [9.93002, -6.46249, 10.5619, 9.37191],
                   [-16.0584, -420.573, -0.472402, -37.2705]])

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

angles = np.arange(6) * 60
rad1 = 7.5 + angles
rad2 = 22.5 + angles
rad3 = 37.5 + angles
rad4 = 52.5 + angles

bend_angle = tuple(zip(rad1, rad2, rad3, rad4))