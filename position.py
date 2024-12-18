# import some libraries
import numpy as np

def Position(x, y, magnet_x, magnet_y, bending_angle):
    # x_diff is the difference between input-x and the magnet's central x
    x_diff = x-magnet_x
    # y_diff is the difference between input-y and the magnet's central y
    y_diff = y-magnet_y
    # rotate coordinate around the Z-axis
    x_particle = (x_diff*np.cos(-np.deg2rad(bending_angle)) -
                  y_diff*np.sin(-np.deg2rad(bending_angle)))
    y_particle = (x_diff*np.sin(-np.deg2rad(bending_angle)) +
                  y_diff*np.cos(-np.deg2rad(bending_angle)))

    coordinate = np.column_stack((x_particle, x_particle/1000, y_particle))  # representing x2, xpos, dy
    return coordinate


# # test the function
# x = 9287.959673
# y = 0.0
# magnet_x = 9232.420563
# magnet_y = 2536.110167
# bending_angle = np.pi / 18  # 45 degrees in radians
# test = Position(x, y, magnet_x, magnet_y, bending_angle)
# # Print the result
# print(test)
