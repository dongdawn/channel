#! /usr/bin/python
#calculate the PMF
import os
import numpy as np                           # The numpy and matplotlib libraries need to be installed!
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import glob
import sys
from scipy.interpolate import spline
font_path = './calibribold.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
leg_path = './Calibri.ttf'
txt_prop = font_manager.FontProperties(fname=leg_path, size=20)
fig = plt.figure(figsize=(6,5)) 
#plt.style.use('ggplot')
sub2 = fig.add_subplot(2,1,1)
rr = np.loadtxt("radius.cs")
sub2.plot(rr[:,0], rr[:,1], color='r', lw=2,alpha=0.7)
p = plt.axvspan(46.8-45.5, 50.9-45.5, facecolor='magenta', alpha=0.072,edgecolor='white')
p = plt.axvspan(38.9-45.5, 43.3-45.5, facecolor='#91F69E', alpha=0.3,edgecolor='#91F69E')
sub2.set_xlim(-35,34)
sub2.set_ylim(0.5,5.5)
#leg=plt.legend(loc=1, labelspacing=0.1, prop={'size': 13.0}, scatterpoints=1, markerscale=1, numpoints=1)
#leg.get_frame().set_linewidth(0.0)
#leg.get_frame().set_alpha(0.1)
sub2.set_ylabel(r'radius ($\mathregular{\AA}$)',fontproperties=font_prop)
sub2.grid(alpha=0.2)
sub2.set_xticklabels([])
for label in (sub2.get_xticklabels() + sub2.get_yticklabels()):
    label.set_fontproperties(font_prop)
    label.set_fontsize(16)
sub1 = fig.add_subplot(2,1,2)
allpmf = np.loadtxt("pmfall_errorbar50000.cso")
x_sm = np.array(allpmf[:,0])
y_sm = np.array(allpmf[:,1])*0.6155
#ymin = np.min(y_sm)
y_sm = y_sm + 3.244
sd = np.array(allpmf[:,2])*0.6155
crystalx = [39.75-45.5,48-45.5,58.25-45.5]
crystaly = [0.042,2.188775714,0.245]
sub1.plot(crystalx,crystaly,"o",color='r',label='Crystal')
sub1.fill_between(x_sm-45.5, y_sm - sd,y_sm + sd, alpha=0.2, edgecolor='#1B2ACC', facecolor='#2D2DFF' )
sub1.plot(x_sm-45.5, y_sm, 'k', color='#1B2ACC', lw=2,label='Glycerol')
p = plt.axvspan(46.8-45.5, 50.9-45.5, facecolor='magenta', alpha=0.072,edgecolor='white')
p = plt.axvspan(38.9-45.5, 43.3-45.5, facecolor='#91F69E', alpha=0.3,edgecolor='#91F69E')
plt.text(47.5-45.5, 7.1, r'SF',color='magenta',rotation='vertical',fontproperties=txt_prop)
plt.text(39.8-45.5, 7.1, r'NPA',color='#00C618', rotation='vertical',fontproperties=txt_prop)
sub1.set_ylabel('PMF (kcal/mol)',fontproperties=font_prop)
sub1.set_xlim(-35,34)
sub1.set_ylim(-0.6,3.6)
#leg=plt.legend(loc=1, labelspacing=0.1, prop={'size': 13.0}, scatterpoints=1, markerscale=1, numpoints=1)
#leg.get_frame().set_linewidth(0.0)
#leg.get_frame().set_alpha(0.1)
sub1.grid(alpha=0.2)
for label in (sub1.get_xticklabels() + sub1.get_yticklabels()):
    label.set_fontproperties(font_prop)
    label.set_fontsize(16)
#p = plt.axvspan(5, 10, facecolor='grey', alpha=0.2,edgecolor='grey')
#p = plt.axvspan(80, 85, facecolor='grey', alpha=0.2,edgecolor='grey')

#sub2.set_yticklabels([])
fig.subplots_adjust(wspace=0, hspace=0)
sub1.set_xlabel(r'pore axis ($\mathregular{\AA}$)',fontproperties=font_prop,labelpad=0)
#plt.savefig('pmf_radius.png',dpi=600)
plt.savefig('pmf_radius.png',dpi=600, bbox_inches='tight')
plt.show()
