# import some libraries
import numpy as np
from decimal import Decimal
# from decimal import Decimal, ROUND_HALF_UP

def Position(x, y, z, magnet_x, magnet_y, bending_angle):
    # x_diff is the difference between input-x and the magnet's central x
    x_diff = Decimal(x-magnet_x)
    # y_diff is the difference between input-y and the magnet's central y
    y_diff = Decimal(y-magnet_y)
    # rotate coordinate around the Z-axis
    x_particle = Decimal(x_diff) * Decimal(np.cos(-np.deg2rad(bending_angle))) - \
                 Decimal(y_diff) * Decimal(np.sin(-np.deg2rad(bending_angle)))
    
    y_particle = Decimal(x_diff) * Decimal(np.sin(-np.deg2rad(bending_angle))) + \
                 Decimal(y_diff) * Decimal(np.cos(-np.deg2rad(bending_angle)))

    # # for calculation accuracy
    # x_particle = Decimal(x_particle).quantize(Decimal('1e-11'), ROUND_HALF_UP)
    # y_particle = Decimal(y_particle).quantize(Decimal('1e-11'), ROUND_HALF_UP)

    coordinate = np.column_stack((float(x_particle), float(y_particle),z))
    
    # print(f'{x_diff*np.cos(-np.deg2rad(bending_angle))}')
    # print(f'{y_diff*np.sin(-np.deg2rad(bending_angle))}')
    # print(f'x_particle is {x_particle/1000}')
    # print(f'mag_rad is {np.deg2rad(bending_angle)}')
    # print(f'explicit calc is {bending_angle/180*np.pi}')
    
    return coordinate

# # test 
# x = 9287.959673
# y = 1700.0
# z = 10.0
# magnet_x = 9232.420563
# magnet_y = 2536.110167
# bending_angle = 7.5
# test = Position(x, y, z, magnet_x, magnet_y, bending_angle)
# # Print the result
# print(test)
# print(type(test))