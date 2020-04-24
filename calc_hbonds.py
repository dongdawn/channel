#!/usr/bin/python

import MDAnalysis
import MDAnalysis.analysis.hbonds
#dirname = ['o6','o9','o11','o12','o16','o17','o18','o22','o24','o25','o26','o30','o31','o32','o33','o34']
dirname = ['o6']
trjname = ['fit_a.xtc','fit_b.xtc','fit_c.xtc','fit_d.xtc']
for d in dirname:
    for t in trjname:
        xtc = "/home/wdd/wdd/glpf_plot/"+d+"/"+t
        print xtc
        outfile = d+"_"+t.split('.')[0]+"_hbonds_as_time.cs"
        wf = open(outfile,'w')
        u = MDAnalysis.Universe('chainA_GOL.pdb', xtc)
        h = MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u, 'resname GOL', distance=3.5, angle=150.0,donors=['OH', 'OG', 'NE2', 'OG1', 'NE', 'N', 'ND1', 'NZ', 'NH1', 'NH2', 'OW', 'ND2', 'SG', 'OH2', 'NE1','O1','O2','O3'],acceptors=['OH', 'OG', 'OD1', 'OD2', 'OG1', 'O', 'ND1', 'NE2', 'OE2', 'OW', 'SG', 'OE1', 'OH2', 'SD','O1','O2','O3'],distance_type='heavy')
        h.run()
        h.generate_table()
#print h.table
        for i in range(len(h.table)):
            for j in range(len(h.table[i])):
                wf.write(str(h.table[i][j])+"       ")
            wf.write("\n")
        wf.close()
