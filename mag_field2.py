import numpy as np
from scipy.spatial import cKDTree

import map
import variables as vs
import position as pos

notrim_map, trim_map = map.map_making()

# make KDTree. x, y, z from trim_map
trim_tree = cKDTree(trim_map[:, :3], compact_nodes=False, balanced_tree=False)
notrim_tree = cKDTree(notrim_map[:, :3], compact_nodes=False, balanced_tree=False)

def mag_field(r):
    x, y, z = r
    pos_mm = np.array([x * 1000, y * 1000, z * 1000])  # convert to mm
    
    # find the nearest magnet
    mag_positions = np.array([vs.magnet_pos_x, vs.magnet_pos_y, vs.magnet_pos_z]).transpose(1, 2, 0)
    mag_distances = np.linalg.norm((mag_positions - pos_mm)**2, axis=2)
    # index of the nearest magnet
    near_mag = np.unravel_index(np.argmin(mag_distances), mag_distances.shape)
    
    # calculate particle position
    particle_pos = pos.Position(
        x * 1000, y * 1000, z * 1000,
        vs.magnet_pos_x[near_mag[0]][near_mag[1]],
        vs.magnet_pos_y[near_mag[0]][near_mag[1]],
        vs.bend_angle[near_mag[0]][near_mag[1]]
    )
    # print(f'particle_pos is {particle_pos}')
    
    # if the nearest magnet is too far
    if abs(particle_pos[0][0]) > 200 or\
       abs(particle_pos[0][1]) > 750 or\
       abs(particle_pos[0][2]) > 18 :
        #    print('too far')
        return np.array([0, 0, 0])

    # note that the order is changed due to the difference between
    # the coordinate of the beam and that of OPERA simulation
    # also unit of position should converted from mm to cm
    else:
        query_point = np.array(
            [particle_pos[0][1]/10, 
             particle_pos[0][0]/10, 
             particle_pos[0][2]/10]
            )
        
        if near_mag[1] in {0, 3}:  # trim magnet
            # find 2 nearest points
            distances, nearest = trim_tree.query(query_point, k=2)
            weights = 1 / (distances + 1e-15)  # weight for each nearest point
            weights /= weights.sum()  # normalized
            B_interp = (weights[0] * trim_map[nearest[0]][[4, 3, 5]] +
                        weights[1] * trim_map[nearest[1]][[4, 3, 5]]) / 10000
            # print(f'particle_pos is {particle_pos}')
            # print(f'weights is {weights}')
            # print(f'nearest point is {trim_map[nearest[0]][[1, 0, 2]]}')
            # print('and')
            # print(f'nearest point is {trim_map[nearest[1]][[1, 0, 2]]}')
            return vs.current_ratio * B_interp

        elif near_mag[1] in {1, 2}:  # no trim magnet
            # find 2 nearest points
            distances, nearest = notrim_tree.query(query_point, k=2)
            weights = 1 / (distances + 1e-15)  # weight for each nearest point
            weights /= weights.sum()  # normalized
            B_interp = (weights[0] * notrim_map[nearest[0]][[4, 3, 5]] +
                        weights[1] * notrim_map[nearest[1]][[4, 3, 5]]) / 10000
            return vs.current_ratio * B_interp

# # test
# x0 = 9.27489285*1000
# y0 = 1850
# z0 = -0.2
# r0 = np.array([x0/1000.0, y0/1000.0, z0/1000.0])
# ans = mag_field(r0)
# print(f'returned magnetic field is {ans}')