#! /usr/bin/python
#get  structure of glycerol along  z scale
from numpy import *                             # The numpy and matplotlib libraries need to be installed!
from MDAnalysis import *        # as well as MDAnalysis
#from matplotlib.pyplot import *
import glob
import sys
import os
curdir  = os.getcwd()
os.chdir("..")
pdb_name = "chainA_GOL.pdb"

#z_begin = 46.8
#out_name = "gol_47_b.pdb"
zmin = 5
zmax = 85
dz = 0.5

xmin = 61                    #scale regions
xmax = 76
ymin = 35.15
ymax = 50.65

s = Universe(pdb_name)
gol = s.select_atoms('resname GOL')
pro = s.select_atoms("protein")

#if the center of mass of gol in the box and return z
def if_contain_gol(g):            #g is the resid of gol 
    sel = s.select_atoms("resid %s" %g)
    xyz = sel.center_of_mass()
    if xyz[0] > xmin and xyz[0] < xmax and xyz[1] > ymin and xyz[1] < ymax and xyz[2] > zmin and xyz[2] < zmax :
        return xyz[2]
    else :
        return -1
W_str = []
for i in arange(zmin,zmax,dz):
    W_str.append('writename'+str(i))
def save_gol_structure(trjname,outname):                #save gol structure in [z,z+dz],
    s.load_new(trjname)
    print "Load  "+trjname+"\n"
    f_log_name = outname + ".log"
    f_log = open(curdir+'/'+f_log_name,'w')
    for ts in s.trajectory:
        sys.stdout.write("Frame  "+str(s.trajectory.frame)+"     Total Frame "+str(len(s.trajectory)))
        sys.stdout.write("\r")
        sys.stdout.flush()
        for gol_index in gol.residues.resids:
            if if_contain_gol(gol_index)!=-1:
                golzz = if_contain_gol(gol_index)
                golzout = int(golzz)
                #if golzz - int(golzz) < 0.5:
                #    golzout = int(golzz)
                #else:
                #    golzout = int(golzz) + 0.5
                pdboutname = outname + "_gol" + str(golzout) + "_" + str(s.trajectory.frame) +".pdb"
                W_str[int(if_contain_gol(gol_index))-zmin] = Writer(curdir+'/'+pdboutname, multiframe = True)
                f_log.write("Frame  "+str(s.trajectory.frame)+"\n")
                f_log.write("Residue Of Glycerol:  "+str(gol_index)+"   in your scale"+"\n")
                gol_str = s.select_atoms('resid %s' %gol_index)
                W_str[int(if_contain_gol(gol_index))-zmin].write(gol_str)
                W_str[int(if_contain_gol(gol_index))-zmin].close()
#    for i in len(W_str):
#        W_str[i].close()
    f_log.close()

#dirindex = [6,11,12,16,22,24,25,26]
dirindex = [6]
for i in dirindex:
    dirname = 'o' + str(i)
    print dirname
    os.chdir("%s" %dirname)
    files = glob.glob("fit*.xtc")
    for f in files:
        out_name = dirname + "/" + dirname + "_" + f.split('.')[0]
        save_gol_structure(f,out_name)
    os.chdir("..")
#os.chdir(curdir)
print "\n",
print "done"
    








