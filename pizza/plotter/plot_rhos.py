#!/usr/bin/env python

import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf


def xilin(x):
      y = 0.804728
      y -= 0.963419*x
      y -= 0.408048*x*x
      y += 0.056507*x*x*x
      y += 0.096629*x*x*x*x
      y -= 0.049933*x*x*x*x*x
      y -= 0.032321*x*x*x*x*x*x
      return y


rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

f = plt.figure(figsize=(10,10))

ax1 = f.add_subplot(221)
ax2 = f.add_subplot(222)
ax3 = f.add_subplot(223)
ax4 = f.add_subplot(224)
axes = [ax1,ax2,ax3,ax4]

c = ['b','r','g']
lt = ['k-.','k','k--']

mass_range_v = ['12.00-13.00','13.00-14.00','14.00-15.00','12.00-15.00']

for mass_range, ax in zip(mass_range_v,axes):

  ax.set_title(mass_range)

  file_proyectada = 'SDSS/DISTORTION_'+mass_range+'.dat'
  data = np.genfromtxt(file_proyectada)
  x = data[:,0]
  y0 = data[:,1]
  y1 = data[:,2]
  y2 = data[:,3]
  ax.plot(x,y0, '-', color='black')
  ax.plot(x,y1, '-', color='red')
  ax.plot(x,y2, '-', color='blue')
  
  file_proyectada = 'SDSS/DISTORTION_rho_200_'+mass_range+'.dat'
  data = np.genfromtxt(file_proyectada)
  x = data[:,0]
  y0 = data[:,1]
  y1 = data[:,2]
  y2 = data[:,3]
  ax.plot(x,y0, '--', color='black')
  ax.plot(x,y1, '--', color='red')
  ax.plot(x,y2, '--', color='blue')



  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim(0.1,30.0)
  ax.set_ylim(1.0,1.0E4)
  ax.set_ylabel(r'$\xi (r)$')
  ax.set_xlabel(r'$r [h^{-1} Mpc]$')
  
  rcParams['legend.frameon']=False
  ax.legend( (r'$\xi_{a}$',r'$\xi_{b}$',r'$\xi_{c}$'), )
  ax.grid(True)

plt.show()
