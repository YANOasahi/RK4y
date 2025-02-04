import numpy as np
import pandas as pd
import glob

# read map files
file_list = glob.glob("./trim_map/*.txt")

# prepare pandas for the data from the files
map_tmp = [pd.read_table(file, header=None) for file in file_list]

# conect pandas
buffer = pd.concat(map_tmp, ignore_index=True)

# convert pandas to ndarray
map = buffer.to_numpy()

# # check
# print(result_array.shape)

# the position of particle
particle = np.array(
    [np.random.randint(-1000, 1000),
     np.random.randint(-200, 200),
     np.random.randint(-20, 20)]
)

print(f'particle position is {particle}')

# pick the nearest point from the map data
distance = [((map[i][0]-particle[0])**2 + (map[i][1]-particle[1]) ** 2 +
            (map[i][2]-abs(particle[2]))**2) for i in range(16800000)]
min = np.argmin(distance)

# check
print(f'min is {min}')
print(f'distance[{min}] is {distance[min]}')

if particle[2] >= 0:
    print(f'the nearest point is ({map[min][0]},{map[min][1]},{map[min][2]})')
    print(f'B is ({map[min][3]},{map[min][4]},{map[min][5]})')
else:
    print(f'the nearest point is ({map[min][0]},{map[min][1]},{-map[min][2]})')
    print(f'B is ({map[min][3]},{-map[min][4]},{map[min][5]})')
