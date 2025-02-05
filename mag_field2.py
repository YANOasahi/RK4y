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
    distances = np.linalg.norm(mag_positions - pos_mm, axis=2)
    # index of the nearest magnet
    near_mag = np.unravel_index(np.argmin(distances), distances.shape)
    
    # calculate particle position
    particle_pos = pos.Position(
        x * 1000, y * 1000, z,
        vs.magnet_pos_x[near_mag[0]][near_mag[1]],
        vs.magnet_pos_y[near_mag[0]][near_mag[1]],
        vs.bend_angle[near_mag[0]][near_mag[1]]
    )
    
    # if the nearest magnet is too far
    if abs(particle_pos[0][0]) > 200 or\
       abs(particle_pos[0][1]) > 1000 or\
       abs(particle_pos[0][2]*1000) > 18 :
        #    print('too far')
           return np.array([0, 0, 0])

    # z>= 0.0 case
    # note that the order is changed due to the difference between
    # the coordinate of the beam and that of OPERA simulation
    elif particle_pos[0][2]>=0:
        query_point = np.array(
            [particle_pos[0][1], 
             particle_pos[0][0], 
             particle_pos[0][2]*1000]
            )
    
        if near_mag[1] in {0, 3}:  # trim magnet
            _, nearest = trim_tree.query(query_point)
            # note that the order is changed due to the difference between
            # the coordinate of the beam and that of OPERA simulation
            return np.array(
                [vs.current_ratio*trim_map[nearest][4], # Bx
                 vs.current_ratio*trim_map[nearest][3], # By
                 vs.current_ratio*trim_map[nearest][5]] # Bz
                )
    
        elif near_mag[1] in {1, 2}: # no trim magnet
            _, nearest = notrim_tree.query(query_point)
            # note that the order is changed due to the difference between
            # the coordinate of the beam and that of OPERA simulation
            return np.array(
                [vs.current_ratio*notrim_map[nearest][4], # Bx 
                 vs.current_ratio*notrim_map[nearest][3], # By 
                 vs.current_ratio*notrim_map[nearest][5]] # Bz
                )
    # z<0.0 case
    # utilize magnetic field symmetry
    elif particle_pos[0][2] < 0:
        query_point = np.array(
            [particle_pos[0][1], 
             particle_pos[0][0], 
             -particle_pos[0][2]*1000]
            )
        # print(f'query point is {query_point}')
    
        if near_mag[1] in {0, 3}:  # trim magnet
            _, nearest = trim_tree.query(query_point)
            # note that the order is changed due to the difference between
            # the coordinate of the beam and that of OPERA simulation
            return np.array(
                [-vs.current_ratio*trim_map[nearest][4], # Bx
                 -vs.current_ratio*trim_map[nearest][3], # By
                 vs.current_ratio*trim_map[nearest][5]]  # Bx
                )
    
        elif near_mag[1] in {1, 2}: # no trim magnet
            _, nearest = notrim_tree.query(query_point)
            # note that the order is changed due to the difference between
            # the coordinate of the beam and that of OPERA simulation
            return np.array(
                [-vs.current_ratio*notrim_map[nearest][4], # Bx
                 -vs.current_ratio*notrim_map[nearest][3], # By
                 vs.current_ratio*notrim_map[nearest][5]]  # Bz
                )
            
# # test
# x0 = 9.27489285*1000
# y0 = 5*1000
# z0 = -0.2
# r0 = np.array([x0/1000.0, y0/1000.0, z0/1000.0])
# ans = mag_field(r0)
# print(f'returned magnetic field is {ans}')