import pandas as pd
import glob

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