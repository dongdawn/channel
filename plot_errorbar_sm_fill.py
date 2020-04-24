#! /usr/bin/python
#calculate the PMF
import os
import numpy as np                           # The numpy and matplotlib libraries need to be installed!
import matplotlib.pyplot as plt
import glob
import sys
from scipy.interpolate import spline
from matplotlib.pyplot import cm
import matplotlib.colors as colors
import matplotlib.font_manager as font_manager
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
cmap = plt.get_cmap('jet_r')
new_cmap = truncate_colormap(cmap, 0.1, 0.9)
font_path = './calibribold.ttf'
leg_path = './Calibri.TTF'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
leg_prop = font_manager.FontProperties(fname=leg_path, size=15.5)
lab_prop = font_manager.FontProperties(fname=leg_path, size=20)
color=new_cmap(np.linspace(0,1,5))
fig = plt.figure(figsize=(6,4))
sub = fig.add_subplot(1,1,1)
zero = [-2.28,-2.71,-2.96,-3.14,-3.27]
for i,c in zip(range(5),color):
    name = "pmfall_errorbar"+str(i+1)+"0000.cso"
    lab = str(i*4.5+4.5)+" $\mathregular{\mu}$s"
    allpmf=np.loadtxt(name)
    x_sm = np.array(allpmf[:,0])-45.5
    y_sm = np.array(allpmf[:,1])*0.6155-zero[i]
    if i == 4:
        c= '#1B2ACC'
    #y_sm = y_sm + 3.244
    #sd = np.array(allpmf[:,2])*0.6155
#x_smooth = np.linspace(x_sm.min(), x_sm.max(), 359)
#y_smooth = spline(x_sm, y_sm, x_smooth)
#plt.fill_between(x_sm, y_sm - sd,y_sm + sd, alpha=0.2, edgecolor='#1B2ACC', facecolor='#2D2DFF' )
    plt.plot(x_sm, y_sm, color=c, lw=2,label = lab)
    #leg=plt.legend(loc=1, labelspacing=0.1, prop={'size': 12.0}, scatterpoints=1, markerscale=1, numpoints=1)
leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
leg.get_frame().set_linewidth(0.0)
leg.get_frame().set_alpha(0.1)

plt.ylabel(r'PMF (kcal/mol)',fontproperties=font_prop)
#plt.xlabel(r'Z',fontsize=15)
p = plt.axvspan(46.8-45.5, 50.9-45.5, facecolor='magenta', alpha=0.072,edgecolor='white')
p = plt.axvspan(38.9-45.5, 43.3-45.5, facecolor='#91F69E', alpha=0.3,edgecolor='white')
plt.text(47.5-45.5, 3.8, r'SF',color='magenta',rotation='vertical',fontproperties=lab_prop)
plt.text(39.8-45.5, 3.8, r'NPA',color='#00C618',rotation='vertical',fontproperties=lab_prop)
#p = plt.axvspan(-35, -30, facecolor='grey', alpha=0.2,edgecolor='grey')
#p = plt.axvspan(29, 34, facecolor='grey', alpha=0.2,edgecolor='grey')
#plt.text(-34, 2.42, r'cytoplasm',rotation='vertical',fontproperties=font_prop)
#plt.text(30, 2.42, r'periplasm',rotation='vertical',fontproperties=font_prop)
#l = plt.axvline(x=46.8, linewidth=4.1, color='hotpink',alpha=0.3)
plt.xlim(-35,34)
plt.ylim(-0.6,4.2)
plt.xlabel(r"pore axis ($\mathregular{\AA}$)",fontproperties=font_prop,labelpad=0)
plt.grid(alpha=0.2)
for label in (sub.get_xticklabels() + sub.get_yticklabels()):
    label.set_fontproperties(font_prop)
    label.set_fontsize(16)
plt.savefig('S1.jpg',dpi=600,bbox_inches='tight')
plt.show()
