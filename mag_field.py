import numpy as np

import map
import variables as vs
import position as pos

notrim_map, trim_map = map.map_making()

def mag_field(r):
    # find the nearest magnet
    x, y, z = r
    distance_mag=[[((vs.magnet_pos_x[i][j]-(x*1000))**2 + 
                    (vs.magnet_pos_y[i][j]-(y*1000))**2 +
                    (vs.magnet_pos_z[i][j]-z*1000)**2) for j in range(4)] for i in range(6)]
    near_mag = np.unravel_index(np.argmin(distance_mag), (6, 4))
    min_distance_mag=distance_mag[near_mag[0]][near_mag[1]]

    # if the nearest magnet is too far,
    # magnetic field will be np.array([0, 0, 0])
    if min_distance_mag > 1000**2 + 200**2 + 18**2:
        return np.array([0, 0, 0])
    
    elif near_mag[1]==0 or near_mag[1]==3:
        # convert particle position to the nearest magnet's coordinate
        particle_pos = pos.Position(x*1000, y*1000, z,
                               vs.magnet_pos_x[near_mag[0]][near_mag[1]],
                               vs.magnet_pos_y[near_mag[0]][near_mag[1]],
                               vs.bend_angle[near_mag[0]][near_mag[1]])
        # find the nearest map data
        pos_in_map = [((notrim_map[i][0]-particle_pos[0][1])**2 + 
                       (notrim_map[i][1]-particle_pos[0][0])**2 +
                       (notrim_map[i][2]-particle_pos[0][2])**2) 
                       for i in range(16800000)]
        nearest = np.argmin(pos_in_map)
        return np.array(
            [notrim_map[nearest][4], 
             notrim_map[nearest][3], 
             notrim_map[nearest][5]]
            )
    
    elif near_mag[1]==1 or near_mag[1]==2:
        # convert particle position to the nearest magnet's coordinate
        particle_pos = pos.Position(x*1000, y*1000, z,
                               vs.magnet_pos_x[near_mag[0]][near_mag[1]],
                               vs.magnet_pos_y[near_mag[0]][near_mag[1]],
                               vs.bend_angle[near_mag[0]][near_mag[1]])
        # find the nearest map data
        pos_in_map = [((trim_map[i][0]-particle_pos[0][1])**2 + 
                       (trim_map[i][1]-particle_pos[0][0])**2 +
                       (trim_map[i][2]-particle_pos[0][2])**2) 
                       for i in range(16800000)]
        nearest = np.argmin(pos_in_map)
        return np.array(
            [trim_map[nearest][4], 
             trim_map[nearest][3], 
             trim_map[nearest][5]]
            )