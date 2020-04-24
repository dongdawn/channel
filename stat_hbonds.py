#!/usr/bin/python
import numpy as np 
from MDAnalysis import * 
pdb_name = "chainA_GOL.pdb"
dirname = ['o6']
#trjname = ['fit_a.xtc','fit_b.xtc','fit_c.xtc','fit_d.xtc']
trjname = ["fit_a.xtc"]
s = Universe(pdb_name)
zmin = 0
zmax = 90
dz = 0.5

cx = 70                 #scale regions
cy = 43
radii = 17
def if_contain_gol(g):            #g is the resid of gol
    sel = s.select_atoms("resid %s" %g)
    xyz = sel.center_of_mass()
    if (xyz[0]-cx)**2 + (xyz[1]-cy)**2 < radii**2 and xyz[2] > zmin and xyz[2] < zmax :
        return xyz[2]
    else :
        return -1
for d in dirname:
    for t in trjname:
        xtc = "/home/wdd/wdd/glpf_plot/"+d+"/"+t
        s.load_new(xtc)
        fname = d+"_"+t.split('.')[0]+"_hbonds_as_time.cs"
        rf = open(fname,'r')
        outname = d+"_"+t.split('.')[0]+"_hbondsz.cs"
        outname2 = d+"_"+t.split('.')[0]+"_hbonds_zonly.cs"
        wf = open(outname,'w')
        wf2 = open(outname2,'w')
        for line in rf: 
            col = line.split()
            time = float(col[0])
            donor_resid = int(col[4])
            acceptor_resid = int(col[7])
            #print int(float(time))
            structure =  s.trajectory[int(time/10)]
            wf.write(str(structure.time)+"     "+col[3]+"     "+col[4]+"     "+col[5]+"     "+col[6]+"     "+col[7]+"     "+col[8])
    
            if donor_resid >= 500 and donor_resid <= 511:
                wf.write("     "+str(if_contain_gol(donor_resid))+"     ")
                wf2.write(str(if_contain_gol(donor_resid))+"     ")
            else:
                wf.write("      0"+"     ")
                wf2.write("0"+"     ")
            if acceptor_resid >= 500 and donor_resid <= 511:
                wf.write("     "+str(if_contain_gol(acceptor_resid))+"\n")
                wf2.write("     "+str(if_contain_gol(acceptor_resid))+"\n")
            else:
                wf.write("      0"+"\n")
                wf2.write("      0"+"\n")
    #for ts in s.trajectory:
    #    if str(s.trajectory.time) == str(time):
    #        print time
    #print s.trajectory.time
        rf.close()
        wf.close()
        wf2.close()
        #countfile = d+"_"+t.split('.')[0]+"_histcount.cs"
        #wf = open(countfile,'w')
        #outname2 = d+"_"+t.split('.')[0]+"_hbonds_zonly.cs"
        #data = np.loadtxt(outname2)
        #donor_data = data[:,0]
        #acceptor_data = data[:,1]
        #bins = np.arange(10,80,0.5)
        #hist_d = np.histogram(donor_data, bins=bins)
        #hist_a = np.histogram(acceptor_data, bins=bins)
        #wf.write("#z          donor          acceptor          total"+"\n")
        #for j in range(len(bins)-1):
        #    wf.write(str(hist_d[1][j]) + "       " + str(hist_d[0][j]) + "       " + str(hist_a[0][j]) + "       " + str(hist_d[0][j]+hist_a[0][j]) + "\n")
        #wf.close()
