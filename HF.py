import numpy as np
import scipy
from scipy.special import erf
from numpy import *
import pandas as pd

def extract_xyz(data_file):

    df = pd.read_csv(data_file, skiprows = 2,
     delim_whitespace= True, names=['atom', 'x', 'y', 'z'])

    atom_type = [i for i in df["atom"]]
    atom_coordinates = [(float(i),float(j),float(k)) for i,j,k in zip(df["x"],df["y"],df["z"])]




    return atom_type, atom_coordinates


print(extract_xyz('HeH.xyz'))
