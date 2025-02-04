import numpy as np
from scipy.spatial import cKDTree

import map
import variables as vs
import position as pos

notrim_map, trim_map = map.map_making()

# make KDTree
trim_tree = cKDTree(trim_map[:, :3])      # x, y, z from trim_map
notrim_tree = cKDTree(notrim_map[:, :3])  # x, y, z from notrim_map

def mag_field(r):
    x, y, z = r
    pos_mm = np.array([x * 1000, y * 1000, z * 1000])  # convert to mm
    
    # find the nearest magnet
    mag_positions = np.array([vs.magnet_pos_x, vs.magnet_pos_y, vs.magnet_pos_z]).transpose(1, 2, 0)
    distances = np.linalg.norm(mag_positions - pos_mm, axis=2)
    near_mag = np.unravel_index(np.argmin(distances), distances.shape)
    min_distance_mag = distances[near_mag]

    # if the nearest magnet is too far
    if min_distance_mag > (1000**2 + 200**2 + 18**2):
        # print('too far!')
        return np.array([0, 0, 0])
    
    # calculate particle position
    particle_pos = pos.Position(
        x * 1000, y * 1000, z,
        vs.magnet_pos_x[near_mag[0]][near_mag[1]],
        vs.magnet_pos_y[near_mag[0]][near_mag[1]],
        vs.bend_angle[near_mag[0]][near_mag[1]]
    )

    # converted particle position to the nearest magnet's coordinate
    query_point = np.array([particle_pos[0][1], particle_pos[0][0], particle_pos[0][2]])
    
    if near_mag[1] in {0, 3}:  # trim magnet
        _, nearest = trim_tree.query(query_point)
        # print('trim magnet')
        return np.array([notrim_map[nearest][4], notrim_map[nearest][3], notrim_map[nearest][5]])
    
    else:  # no trim magnet
        _, nearest = notrim_tree.query(query_point)
        # print('no trim magnet')
        return np.array([trim_map[nearest][4], trim_map[nearest][3], trim_map[nearest][5]])
