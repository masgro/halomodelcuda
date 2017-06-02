#!/usr/bin/env python

import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf,gamma
import os.path
from scipy import integrate
import functions as funcs

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

f = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(4,1,)
ax = f.add_subplot(gs[0:3,0])
ax1 = f.add_subplot(gs[3,0],sharex=ax)

gs.update(hspace=0.0)
plt.setp(ax.get_xticklabels(), visible=False)

folder = '../outputs/chi/2ndpass/'
mag = 17
mass_range = '13.33-13.74'
bc = 0.89
ab = 0.50   
al_b = 0.80
al_c = 1.30

funcs.plot_halo_model(folder,mag,mass_range,bc,ab,al_b,al_c,ax,ax1)

#bc = np.linspace(0.7,0.9,5)   
#ab = np.linspace(0.7,0.9,5)   
#
#bcv = np.repeat(bc,5)
#abv = np.tile(ab,5)
#
al_b = np.linspace(0.3,1.3,21)
al_c = np.linspace(0.7,1.7,21)

al_bv = np.repeat(al_b[5],21)
al_cv = np.tile(al_c,21)

for i in range(21):
    #funcs.plot_halo_model(folder,mag,mass_range,bcv[i],abv[i],al_bv[i],al_cv[i],ax,ax1)
    funcs.plot_halo_model(folder,mag,mass_range,0.89,0.50,al_bv[i],al_cv[i],ax,ax1)

filename = 'facu/mockscorpio/xiproy_-'+'%2d'%mag+'_'+mass_range+'_20bin_100mp.dat'
if os.path.isfile(filename):
    r,xi,xiper,xipar=np.genfromtxt(filename,unpack=True,usecols=(0,1,2,3))
    ax.plot(r,xi,'-ok',lw=1,label=r"Proyectada desde $\xi(\sigma,\pi)$")
    ax.plot(r,xiper,'-ob',lw=1,label=r"Proyectada $\perp$")
    ax.plot(r,xipar,'-or',lw=1,label=r"Proyectada $\parallel$")

    ax1.plot(r,xi/xi,'-ok')
    ax1.plot(r,xiper/xi,'-ob',lw=1)
    ax1.plot(r,xipar/xi,'-or',lw=1)


ax1.set_xscale('log')
ax1.set_xlabel(r'$r [h^{-1} Mpc]$',fontsize=18)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,100.0)
ax.set_ylim(1.0,5.0E3)
ax.set_ylabel(r'$\omega (r_{p})$',fontsize=18)

ax.set_title('Mass range: %s Magnitude: %2d'%(mass_range,mag))
rcParams['legend.frameon']=True
ax.legend()
ax.grid(True)

plt.show()
