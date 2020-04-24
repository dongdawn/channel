#!/usr/bin/env python
import os
from numpy import *                             # The numpy and matplotlib libraries need to be installed!
from MDAnalysis import *        # as well as MDAnalysis
from matplotlib.pyplot import *
import sys

#dirname = ['o6','o9','o11','o12','o16','o17','o18','o22','o24','o25','o26','o30','o31','o32','o33','o34']
#dirname = ['o6','o9','o11','o12','o16','o18','o22','o24','o26','o30','o31','o32','o33','o34','v1','v3','v4','v6','v7','v8','v9','v10','v13','v14','v16','v17','v19','v20','v23','v24','v25','v2','v5','v11','v12','v15','v18','v27','v28','v21','v26','v29','v30']
dirname = ['v1','v3','v4','v6','v7','v8','v9','v10','v13','v14','v16','v17','v19','v20','v23','v24','v25','v2','v5','v11','v12','v15','v18','v27','v28','v21','v26','v29','v30']
#dirname = ['o25']
trjname = ['fit_a.xtc','fit_b.xtc','fit_c.xtc','fit_d.xtc']
#trjname = ['fit_a.xtc','fit_c.xtc','fit_d.xtc']
pdb_name = "chainA_GOL.pdb"
s = Universe(pdb_name)
gol = s.select_atoms('resname GOL')
#gol = s.select_atoms('resname SOL and cyzone 18 44 -44 protein')
pro = s.select_atoms("protein")

cx = 70                 #scale regions
cy = 43
radii = 17
cz = 47.25
def if_contain_gol(g):            #g is the resid of gol 
    sel = s.select_atoms("resid %s" %g)
    xyz = sel.center_of_mass()
    if (xyz[0]-cx)**2 + (xyz[1]-cy)**2 < radii**2 and xyz[2] > cz and xyz[2] < cz+0.5:
        return xyz[2]
    else :
        return 0
def normal(p1,p2,p3):               #return the normal of three points
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    # These two vectors are in the plane
    v1 = p2 - p3
    v2 = p2 - p1
    # the cross product is a vector normal to the plane
    cp = np.cross(v1, v2)
    return cp
def cosVector(x,y):
    if(len(x)!=len(y)):
        print('error input,x and y is not in the same space')
        return;
    x = array(x)
    y = array(y)
    lx = sqrt(x.dot(x))
    ly = sqrt(y.dot(y))
    return (x.dot(y)/(lx*ly))
def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = m.sqrt(XsqPlusYsq + z**2)               # r
    elev = m.atan2(z,m.sqrt(XsqPlusYsq))     # theta
    az = m.atan2(y,x)                           # phi
    return r, elev, az
for dn in dirname:
    for tn in trjname:
        out_name = dn + '_chain' + tn.split('.')[0].split('_')[1] + '_sf_dihed.dat'
        out_name2 = dn + '_chain' + tn.split('.')[0].split('_')[1] + '_sf_tiltangle.dat'
        trj_name = '/home/wdd/glpf/density_calculate_version_c/' + dn + '/' + tn
        print trj_name
        fw2 = open(out_name2,"w")
        fw2.write("#Frame  48theta  48phi  200theta 200phi    z"+"\n")
        fw = open(out_name,"w")
        fw.write("#Frame  48CA-CB-CG-CD2  200CA-CB-CG-CD2  206CA-CB-CG-CD 206CB-CG-CD-NE  206CG-CD-NE-CZ    Z"+"\n")
        s.load_new(trj_name)
        for ts in s.trajectory[::10]:
            sys.stdout.write("Frame  "+str(s.trajectory.frame)+"     Total Frame "+str(len(s.trajectory)))
            sys.stdout.write("\r")
            sys.stdout.flush()
            for gol_index in gol.residues.resids:
                if if_contain_gol(gol_index) != 0:
                    d0 = s.select_atoms('bynum 623','bynum 625',' bynum 628', ' bynum 634')
                    d1 = s.select_atoms('bynum 2907',' bynum 2909',' bynum 2912','bynum 2919')
                    d2 = s.select_atoms('bynum 2992',' bynum 2994',' bynum 2997','bynum 3000')
                    d3 = s.select_atoms('bynum 2994',' bynum 2997',' bynum 3000','bynum 3003')
                    d4 = s.select_atoms('bynum 2997',' bynum 3000',' bynum 3003','bynum 3005')
                    dih0 = d0.dihedral.value()
                    dih1 = d1.dihedral.value()
                    dih2 = d2.dihedral.value()
                    dih3 = d3.dihedral.value()
                    dih4 = d4.dihedral.value()
                    trp1 = s.select_atoms('bynum 635')
                    trp2 = s.select_atoms('bynum 637')
                    trp3 = s.select_atoms('bynum 641')
                    phe1 = s.select_atoms('bynum 2915')
                    phe2 = s.select_atoms('bynum 2917')
                    phe3 = s.select_atoms('bynum 2921')
                    trpp1 = trp1.positions[0]
                    trpp2 = trp2.positions[0]
                    trpp3 = trp3.positions[0]
                    phep1 = phe1.positions[0]
                    phep2 = phe2.positions[0]
                    phep3 = phe3.positions[0]
                    normal_trp = normal(trpp1,trpp2,trpp3)
                    normal_phe = normal(phep1,phep2,phep3)
                    tr,ttheta,tphi = cart2sph(normal_trp[0],normal_trp[1],normal_trp[2])
                    pr,ptheta,pphi = cart2sph(normal_phe[0],normal_phe[1],normal_phe[2])
                    fw.write(str(s.trajectory.frame)+"      "+str(dih0)+"      "+str(dih1)+"      "+str(dih2)+"      "+str(dih3)+"      "+str(dih4)+"      "+str(if_contain_gol(gol_index)))
                    fw.write("\n")
                    fw2.write(str(s.trajectory.frame)+"      "+str(90-ttheta*180/np.pi)+"      " + str(tphi*180/np.pi)+"      "+str(90-ptheta*180/np.pi)+"      "+str(pphi*180/np.pi)+"      "+str(if_contain_gol(gol_index)))
                    fw2.write("\n")
        fw.close()
        fw2.close()
        print "\n"
        print "done"

