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
leg_prop = font_manager.FontProperties(fname=leg_path, size=15)
txt_prop = font_manager.FontProperties(fname=leg_path, size=20)
color=new_cmap(np.linspace(0,1,5))
fig = plt.figure(figsize=(6,4))
sub = fig.add_subplot(1,1,1)
allpmf = np.loadtxt("alltrans_pro1.cso")
x_sm = np.array(allpmf[:,0])
y_sm = np.array(allpmf[:,1])*100
sd = np.array(allpmf[:,2])*100
oridata = np.loadtxt("ave_count_errorbar_orien.cso")
x_o = np.array(oridata[:,0])
y_o = np.array(oridata[:,1])*100
sd_o = np.array(oridata[:,2])*100

condata = np.loadtxt("ave_count_errorbar_confor.cso")
x_c = np.array(condata[:,0])
y_c = np.array(condata[:,1])*100
sd_c = np.array(condata[:,2])*100
#ymin = np.min(y_sm)
sub.fill_between(x_sm-45.5, y_sm - sd,y_sm + sd, alpha=0.2, edgecolor='#1B2ACC', facecolor='#2D2DFF' )
sub.plot(x_sm-45.5, y_sm, 'k', color='#1B2ACC', lw=2,label='all')

sub.fill_between(x_o-45.5, y_o - sd_o,y_o + sd_o, alpha=0.2, edgecolor='red', facecolor='red' )
sub.plot(x_o-45.5, y_o, 'k', color='red', lw=2,label='orientation')

sub.fill_between(x_c-45.5, y_c - sd_c,y_c + sd_c, alpha=0.2, edgecolor='orange', facecolor='orange' )
sub.plot(x_c-45.5, y_c, 'k', color='orange', lw=2,label='conformation')

p = plt.axvspan(46.8-45.5, 50.9-45.5, facecolor='magenta', alpha=0.072,edgecolor='white')
p = plt.axvspan(38.9-45.5, 43.3-45.5, facecolor='#91F69E', alpha=0.3,edgecolor='#91F69E')
p = plt.axvspan(2, 2.25, facecolor='gray', edgecolor='gray',alpha=0.3)
plt.text(47.8-45.5, 45.3, r'SF',color='magenta',rotation='vertical',fontproperties=txt_prop)
plt.text(39.4-45.5, 45.3, r'NPA',color='#00C618', rotation='vertical',fontproperties=txt_prop)
sub.set_ylabel('probability (%)',fontproperties=font_prop)
sub.set_xlim(-35,34)
sub.set_ylim(0,50)
leg=plt.legend(loc=1, labelspacing=0.1, prop=leg_prop, scatterpoints=1, markerscale=1, numpoints=1,handlelength=1.5)
leg.get_frame().set_linewidth(0.0)
leg.get_frame().set_alpha(0.1)
for label in (sub.get_xticklabels() + sub.get_yticklabels()):
    label.set_fontproperties(font_prop)
    label.set_fontsize(16)
#p = plt.axvspan(5, 10, facecolor='grey', alpha=0.2,edgecolor='grey')
#p = plt.axvspan(80, 85, facecolor='grey', alpha=0.2,edgecolor='grey')

#sub2.set_yticklabels([])
fig.subplots_adjust(wspace=0, hspace=0)
sub.set_xlabel(r'pore axis ($\mathregular{\AA}$)',fontproperties=font_prop,labelpad=0)
#plt.savefig('pmf_radius.png',dpi=600)
plt.savefig('alltranspro_all.jpg',dpi=600, bbox_inches='tight')
plt.show()
