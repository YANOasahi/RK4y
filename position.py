# import some libraries
import numpy as np

# import variables defined by myself
import variables_position as vp


def Position(x, y, magnet_x, magnet_y, bending_angle):
    # x_diff is the difference between input-x and the magnet's central x
    x_diff = x-magnet_x
    # y_diff is the difference between input-y and the magnet's central y
    y_diff = y-magnet_y
    # rotate coordinate around the Z-axis
    x_particle = (x_diff*np.cos(bending_angle) -
                  y_diff*np.sin(bending_angle))/1000
    y_particle = (x_diff*np.sin(bending_angle) +
                  y_diff*np.cos(bending_angle))

    coordinate = np.column_stack((x_particle, y_particle))
    return coordinate
    

# # test the function
# x = 1.0
# y = 2.0
# magnet_x = 0.5
# magnet_y = 1.0
# bending_angle = np.pi / 4  # 45 degrees in radians
# test = Position(x, y, magnet_x, magnet_y, bending_angle)
# # Print the result
# print(test)
