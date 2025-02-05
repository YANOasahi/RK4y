import pandas as pd
import glob
import matplotlib.pyplot as plt

def map_making():
    # read map files
    notrim_files = glob.glob("./notrim_map/*.txt")
    trim_files = glob.glob("./trim_map/*.txt")
    # prepare pandas for the data from the files
    notrim_tmp = [pd.read_table(file, header=None) for file in notrim_files]
    trim_tmp = [pd.read_table(file, header=None) for file in trim_files]
    # conect pandas
    notrim_buffer = pd.concat(notrim_tmp, ignore_index=True)
    trim_buffer = pd.concat(trim_tmp, ignore_index=True)
    # convert pandas to ndarray
    notrim = notrim_buffer.to_numpy()
    trim = trim_buffer.to_numpy()
    return notrim,trim

notrim_map, trim_map=map_making()
print(notrim_map)
print(notrim_map.shape)
print('*******************')
print(trim_map)
print(trim_map.shape)
# 3-D scatter plot
fig1_1, ax = plt.subplots(figsize=(9, 9), subplot_kw={'projection': '3d'})
ax.scatter(notrim_map[:,0], notrim_map[:,1], -notrim_map[:,5])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Bz")
ax.set_title("notrim")
fig1_2, bx = plt.subplots(figsize=(9, 9), subplot_kw={'projection': '3d'})
bx.scatter(trim_map[:,0], trim_map[:,1], -trim_map[:,5])
bx.set_xlabel("x")
bx.set_ylabel("y")
bx.set_zlabel("Bz")
bx.set_title("trim")

plt.legend()
plt.show()