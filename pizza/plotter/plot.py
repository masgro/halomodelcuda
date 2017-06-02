#!/usr/bin/env python

import os.path
import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf


import matplotlib.gridspec as gridspec

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

# Three subplots sharing both x/y axes
f, axes = plt.subplots(2,2, sharex=True, sharey=True, figsize=(12,12))
f.subplots_adjust(hspace=0,wspace=0.0)

plt.setp([a.get_xticklabels() for a in f.axes[:-4]], visible=False)
plt.setp([a.get_yticklabels() for a in f.axes[1:3]], visible=False)
plt.setp([a.get_yticklabels() for a in f.axes[5:7]], visible=False)

c = ['b','r','g']
lt = ['k-.','k','k--']

nlines = 50
ndir = 2
ym1h = np.zeros((ndir, nlines))
er1h = np.zeros((ndir, nlines))
ym2h = np.zeros((ndir, nlines))
er2h = np.zeros((ndir, nlines))
ym   = np.zeros((ndir, nlines))
er   = np.zeros((ndir, nlines))

FACTOR = 0.6

mag = 17
mass_range_v = ['12.50-12.97','12.97-13.33','13.33-13.74','13.74-15.00']
kind = 'CG'

for mass_range, ax in zip(mass_range_v,axes.reshape(4)):

  ax.set_title(mass_range)

  filename1h = '../output/funcorr_1h_'+mass_range+'_'+'%2d'%mag+'.'+kind+'.dat'
  filename2h = '../output/funcorr_2h_'+mass_range+'_'+'%2d'%mag+'.'+kind+'.dat'
  
  if os.path.isfile(filename1h) & os.path.isfile(filename2h):
    data = np.genfromtxt(filename1h)
    x = data[:,0]
    ym1h[0,:] = data[:,1]*FACTOR
    er1h[0,:] = data[:,2]
    ym1h[1,:] = data[:,3]*FACTOR
    er1h[1,:] = data[:,4]
  
    data = np.genfromtxt(filename2h)
    x = data[:,0]
    ym2h[0,:] = data[:,1]*FACTOR
    er2h[0,:] = data[:,2]
    ym2h[1,:] = data[:,3]*FACTOR
    er2h[1,:] = data[:,4]
    
    ym = ym1h + ym2h
    er = np.sqrt(er1h*er1h + er2h*er2h)
    
    y = ym[0,:]
    #color = 'blue'
    #ax.plot(x,y, '-', color=color)
    
    y = ym[1,:]
    color = 'black'
    ax.plot(x,y, '-',lw=1.5, color=color,label='Modelo')

    #y = ym1h[0,:]
    #ax.plot(x,y, ':', color='black')
    y = ym1h[1,:]
    ax.plot(x,y, ':', color='black')

    #y = ym2h[0,:]
    #ax.plot(x,y, ':', color='black')
    y = ym2h[1,:]
    ax.plot(x,y, ':', color='black')
  else :
    print 'NO EXISTE EL ARCHIVO: ',filename1h
  

  filename = 'facu/mockscorpio/xiproy_-'+'%2d'%mag+'_'+mass_range+'.dat'
  if os.path.isfile(filename):
      r,xi=np.genfromtxt(filename,unpack=True,usecols=(0,1))
      ax.plot(r,xi,'-oy',lw=3,label='MockScorpio')
  else :
      print 'NO EXISTE EL ARCHIVO: ',filename

  suffix = mass_range+'_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (90,0.99,0.99,1.4,1.4,mag)
  filename = '../../cross/plotter/proyectada_integrada_'+suffix+'.dat'
  if os.path.isfile(filename):
      r,xi,t1h,t2h=np.genfromtxt(filename,unpack=True,usecols=(0,1,2,3))
      ax.plot(r,xi,'-r',lw=2,label="Modelo 3D proyectado")
      ax.plot(r,t1h,':r',lw=2)
      ax.plot(r,t2h,':r',lw=2)
  else :
      print 'No existe file ',filename

  #suffix = mass_range
  #filename = '../../cross/plotter/integrado_xibox.-21.'+suffix+'.dat'
  #if os.path.isfile(filename):
  #    r,xi=np.genfromtxt(filename,unpack=True,usecols=(0,1))
  #    ax.plot(r,xi,'-b',lw=2,label="Integral Modelo 3D")
  #else :
  #    print filename

  #filename = '../../cross/plotter/integrado_xibox.-17.dat'
  #if os.path.isfile(filename):
  #    r,xi=np.genfromtxt(filename,unpack=True,usecols=(0,1))
  #    ax.plot(r,xi,'-b',lw=2,label="Integral 3D")
  #else :
  #    print filename

  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlim(0.1,30.0)
  ax.set_ylim(1.0,1.0E4)
  ax.set_ylabel(r'$\xi (r)$')
  ax.set_xlabel(r'$r [h^{-1} Mpc]$')
  
  rcParams['legend.frameon']=False
  ax.legend()
  ax.grid(True)

plt.show()
