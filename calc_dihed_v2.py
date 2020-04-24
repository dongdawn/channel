#! /usr/bin/python
from MDAnalysis import *
import numpy as np
import os
#for i in np.arange(6,85):
for i in [47.5]:
    outname = 'dihed' + str(i) + '.cs'
    f = open(outname,'w')
    f.write("# O1-C1-C2-O2  O3-C3-C2-O2  OH1-O1-C1-C2  OH2-O2-C2-C1  OH3-O3-C3-C2" + "\n")
    pdb_name = "all_new_gol" + str(i) + ".pdb"
    s = Universe(pdb_name)
    for ts in s.trajectory:
        d0 = s.select_atoms('bynum 4','bynum 1',' bynum 6', ' bynum 8')
        d1 = s.select_atoms('bynum 13',' bynum 10',' bynum 6','bynum 8')
        d2 = s.select_atoms('bynum 5',' bynum 4',' bynum 1','bynum 6')
        d3 = s.select_atoms('bynum 9',' bynum 8',' bynum 6','bynum 1')
        d4 = s.select_atoms('bynum 14',' bynum 13',' bynum 10','bynum 6')
        dih0 = d0.dihedral.value()
        dih1 = d1.dihedral.value()
        dih2 = d2.dihedral.value()
        dih3 = d3.dihedral.value()
        dih4 = d4.dihedral.value()
        f.write(str(dih0) + "      " + str(dih1) + "      " + str(dih2) + "      " +str(dih3) + "      " + str(dih4) + "      ")
        f.write('\n')
    f.close()
