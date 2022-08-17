#!/usr/bin/python

import numpy as np
import diskcache
import os, sys

sys.path.append(os.getcwd())

#read in names of xyz files
readfile = 'names.txt'
with open(readfile, 'r') as f:
    systems = f.readlines()

data = diskcache.Index('/scr/weber/G2_db')

print(systems)

el_to_atom = { 'H' : 1,
        'He' : 2,
        'Li' : 3,
        'Be' : 4,
        'B' : 5,
        'C' : 6,
        'N' : 7,
        'O' : 8,
        'F' : 9,
        'Ne' : 10,
        'Na' : 11,
        'Mg' : 12,
        'Al' : 13,
        'Si' : 14,
        'P' : 15,
        'S' : 16,
        'Cl' : 17,
        'Ar' : 18}

i = 0
for mol in systems:
    mol = mol.strip('\n')
    elements = np.genfromtxt(mol,dtype='str',skip_header=2,usecols=0)
    atomicNumbers = [el_to_atom[el] for el in elements]
    xyz = np.loadtxt(mol, dtype=np.float64, skiprows=2, usecols=(1,2,3))
    data[i] = {'name' : mol, 'positions' : xyz, 'elements' : elements, 'atomicNumbers' : atomicNumbers}
    i += 1
