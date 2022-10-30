import numpy as np
import scipy
from scipy.special import erf
from numpy import *
import pandas as pd


### Extract data from xyz file ###

def extract_xyz(data_file):

    dg = pd.read_csv(data_file, nrows = 1, names = ['num_atom'])
    df = pd.read_csv(data_file, skiprows = 2,
     delim_whitespace= True, names=['atom', 'x', 'y', 'z'])

    num_atoms = dg["num_atom"][0]
    atom_type = [i for i in df["atom"]]
    atom_coordinates = [np.array([float(i),float(j),float(k)]) for i,j,k in zip(df["x"],df["y"],df["z"])]



    return num_atoms, atom_type, atom_coordinates


### Product of gaussians ###

def gauss_product(gauss_A,gauss_B):
    ## Pass the exponent and centre as a tuple
    a, Ra = gauss_A
    b, Rb = gauss_B
    p = a + b
    diff = np.linalg.norm(Ra-Rb)**2
    N = (4*a*b/(pi**2))**0.75
    K = N*exp(-a*b/p*diff)
    Rp = (a*Ra + b*Rb)/p

    return p, diff, K, Rp

## overlap integral
def overlap(A,B):
    p, diff, K, Rp = gauss_product(A,B)
    prefactor = (pi/p)**1.5
    return prefactor*K

## kinetic integral
def kinetic_integral(A,B):
    p, diff, K, Rp = gauss_product(A,B)
    prefactor = (pi/p)**1.5

    a, Ra = A
    b, Rb = B

    reduced_exponent = a*b/p
    return reduced_exponent*(3-2*reduced_exponent*diff)*prefactor*K

## The Boys function
def Fo(t):
    if (t == 0):
        return 1
    else:
        return (0.5*(pi/t)**0.5)*erf(t**0.5)


## Nuclear-electron integral(the C potential)
def potential(A,B, atom_idx):
    p, diff, K, Rp = gauss_product(A,B)
    Rc = atom_coordinates[atom_idx]
    Zc = charge_dict[atoms[atom_idx]]

    return (-2*pi*Zc/p)*K*Fo(p*np.linalg.norm(Rp - Rc)**2)

## Two electron integral
def multi(A,B,C,D):
    p, diff_ab, K_ab, Rp = gauss_product(A,B)
    q, diff_cd, K_cd, Rq = gauss_product(C,D)
    multi_prefactor = 2*pi**2.5*(p*q*(p+q)**0.5)*-1

    return multi_prefactor*K_ab*K_cd*Fo(p*q/(p+q)*np.linalg.norm(Rp-Rq)**2)


number_of_atoms = 0
atom_types = []
atom_coordinates = []

number_of_atoms, atom_types, atom_coordinates = extract_xyz('HeH.xyz')

#print(atom_coordinates[1])
#print(np.linalg.norm(atom_coordinates[1]-atom_coordinates[0]))



### Defining the basis ####

STOng = 3
zeta_dict = {'H':[1.24], 'He':[2.095], 'Li':[2.69]}
max_quantum_number = {'H':1,'He':1}

D = np.array([[0.444635, 0.535328, 0.154329],[0.700115,0.399513,-0.0999672]])

alpha = np.array([[0.109818, 0.405771, 2.22766],[0.0751386, 0.231031, 0.994203]])


B = 0
for atom in atom_types:
    B += max_quantum_number[atom]



#Number of electrons

N = 2


## Charges

charge_dict =   {'H':1, 'He':2}
