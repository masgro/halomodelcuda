#!/usr/bin/env python

import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf
import os.path

def modelo(mass_range,ax):

    escalas = np.zeros((3,2))
    sigmas = np.zeros((3,2))

    nlines = 50
    ndir = 3
    ym1h = np.zeros((ndir, nlines))
    er1h = np.zeros((ndir, nlines))
    ym2h = np.zeros((ndir, nlines))
    er2h = np.zeros((ndir, nlines))
    ym   = np.zeros((ndir, nlines))
    er   = np.zeros((ndir, nlines))

    if(mass_range == '12.50-12.97'):
        sbc = '0.10' 
        sab = '0.20'
        sal_b = '0.60'
        sal_c = '1.00'
        factor_1h = 1.2
        factor_2h = 1.15
        escalas[0,1] = 0.1        ; sigmas[0,1] = 0.96548314
        escalas[0,0] = 0.32641641 ; sigmas[0,0] = 6.99992726
        escalas[1,1] = 0.55045432 ; sigmas[1,1] = 4.02410762
        escalas[1,0] = 0.93740779 ; sigmas[1,0] = 0.74814675
        escalas[2,1] = 0.23862464 ; sigmas[2,1] = 1.73605431
        escalas[2,0] = 0.32167312 ; sigmas[2,0] = 6.58714271

    if(mass_range == '12.97-13.33'):
        sbc = '0.25'
        sab = '0.10'
        sal_b = '0.90'
        sal_c = '0.70'
        factor_1h = 1.2
        factor_2h = 1.15
        escalas[0,1] = 0.53343258; sigmas[0,1] = 3.29455895
        escalas[0,0] = 0.41179681; sigmas[0,0] = 3.47114219
        escalas[1,1] = 0.44917009; sigmas[1,1] = 3.12758295
        escalas[1,0] = 0.35662938; sigmas[1,0] = 4.23206749
        escalas[2,1] = 0.47374043; sigmas[2,1] = 2.25528555
        escalas[2,0] = 0.4137295 ; sigmas[2,0] = 4.0970165 

    if(mass_range == '13.33-13.74'):
        sbc = '0.35' 
        sab = '0.20'
        sal_b = '0.70'
        sal_c = '0.30'
        factor_1h = 1.2
        factor_2h = 1.0

        escalas[0,0] = 0.66
        escalas[1,0] = 0.89
        escalas[2,0] = 0.76
        escalas[0,1] = 0.50
        escalas[1,1] = 1.12
        escalas[2,1] = 1.067

        sigmas[0,0] =  5.87
        sigmas[1,0] =  1.56
        sigmas[2,0] =  3.87
        sigmas[0,1] =  1.558
        sigmas[1,1] =  3.46
        sigmas[2,1] =  4.511

    if(mass_range == '13.74-15.00'):
        sbc = '0.50' 
        sab = '0.20'
        sal_b = '0.55'
        sal_c = '0.30'
        factor_1h = 1.2
        factor_2h = 1.0

        escalas[0,0] = 1.22901972
        escalas[1,0] = 1.42679342
        escalas[2,0] = 1.23729639
        escalas[0,1] = 0.06709902
        escalas[1,1] = 1.69025408 
        escalas[2,1] = 0.77660862
        sigmas[0,0] = 4.2920803
        sigmas[1,0] = 3.8874773
        sigmas[2,0] = 4.38157754
        sigmas[0,1] = 0.65937736
        sigmas[1,1] = 3.47046722
        sigmas[2,1] = 1.88762472

    angulo = '45'
    kind = 'CG'
    mag = 17
    folder = '../outputs/'
    suffix = mass_range+'_'+angulo+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+'%2d'%mag+'.'+kind+'.dat'
    
    filename1h = folder+'funcorr_1h_'+suffix
    filename2h = folder+'funcorr_2h_'+suffix
    
    if os.path.isfile(filename1h):
        data = np.genfromtxt(filename1h,unpack=True)
        x = data[0]
        ym1h[0,:] = data[1]*factor_1h
        er1h[0,:] = data[2]*factor_1h
        ym1h[1,:] = data[3]*factor_1h
        er1h[1,:] = data[4]*factor_1h
        ym1h[2,:] = data[5]*factor_1h
        er1h[2,:] = data[6]*factor_1h
    
    if os.path.isfile(filename2h):
        data = np.genfromtxt(filename2h,unpack=True)
        x = data[0]
        ym2h[0,:] = data[1]*factor_2h
        er2h[0,:] = data[2]*factor_2h
        ym2h[1,:] = data[3]*factor_2h
        er2h[1,:] = data[4]*factor_2h
        ym2h[2,:] = data[5]*factor_2h
        er2h[2,:] = data[6]*factor_2h
    
    
    for i in range(3):
        ym1h[i]*=(erf((np.log10(escalas[i,0])-np.log10(x))*sigmas[i,0])*0.5+0.5)
        ym2h[i]*=(erf((np.log10(x)-np.log10(escalas[i,1]))*sigmas[i,1])*0.5+0.5)
    
    
    ym = ym1h + ym2h
    er = np.sqrt(er1h*er1h + er2h*er2h)
    
    y = ym[0,:]
    error = er[0,:]
    color = 'red'
    ax.plot(x,y, '-', lw=2,color=color)
    
    y = ym[1,:]
    error = er[1,:]
    color = 'blue'
    ax.plot(x,y, '-', lw=2,color=color)
    
    y = ym[2,:]
    error = er[2,:]
    color = 'black'
    ax.plot(x,y,'-',lw=2,color=color,label="Modelo")
    
    #cocientes
    #interp = interp1d(x,ym[0,:])
    #color = 'red'
    #ax1.plot(r,interp(r)/xi, '-', lw=2,color=color)
    #color = 'blue'
    #interp = interp1d(x,ym[1,:])
    #ax1.plot(r,interp(r)/xi, '-', lw=2,color=color)
    #color = 'black'
    #interp = interp1d(x,ym[2,:])
    #ax1.plot(r,interp(r)/xi, '-', lw=2,color=color)
    #############################################################################################



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
gs = gridspec.GridSpec(8,2)
ax11 = f.add_subplot(gs[0:3,0])
ax12 = f.add_subplot(gs[3,0],sharex=ax11)
ax21 = f.add_subplot(gs[4:7,0])
ax22 = f.add_subplot(gs[7,0],sharex=ax21)
ax31 = f.add_subplot(gs[0:3,1],sharey=ax11)
ax32 = f.add_subplot(gs[3,1],sharex=ax31,sharey=ax12)
ax41 = f.add_subplot(gs[4:7,1],sharey=ax21)
ax42 = f.add_subplot(gs[7,1],sharex=ax41,sharey=ax22)

gs.update(hspace=0.0)
gs.update(wspace=0.0)

plt.setp(ax11.get_xticklabels(), visible=False)
plt.setp(ax21.get_xticklabels(), visible=False)
plt.setp(ax31.get_xticklabels(), visible=False)
plt.setp(ax41.get_xticklabels(), visible=False)
plt.setp(ax12.get_xticklabels(), visible=False)
plt.setp(ax32.get_xticklabels(), visible=False)

plt.setp(ax31.get_yticklabels(), visible=False)
plt.setp(ax32.get_yticklabels(), visible=False)
plt.setp(ax41.get_yticklabels(), visible=False)
plt.setp(ax42.get_yticklabels(), visible=False)

ax11.set_xscale('log')
ax21.set_xscale('log')
ax31.set_xscale('log')
ax41.set_xscale('log')
ax12.set_xscale('log')
ax22.set_xscale('log')
ax32.set_xscale('log')
ax42.set_xscale('log')
ax11.set_yscale('log')
ax21.set_yscale('log')
ax31.set_yscale('log')
ax41.set_yscale('log')

ax11.set_ylabel(r'$\omega (r_{p})$',fontsize=18)
ax21.set_ylabel(r'$\omega (r_{p})$',fontsize=18)
ax22.set_xlabel(r'$r [h^{-1} Mpc]$',fontsize=18)
ax42.set_xlabel(r'$r [h^{-1} Mpc]$',fontsize=18)

ax12.set_ylim(0.0,2.0)
ax22.set_ylim(0.0,2.0)
ax32.set_ylim(0.0,2.0)
ax42.set_ylim(0.0,2.0)

ax11.set_xlim(0.1,23.0)
ax11.set_ylim(6.0,5.0E3)
ax21.set_xlim(0.1,23.0)
ax21.set_ylim(6.0,5.0E3)
ax31.set_xlim(0.1,23.0)
ax31.set_ylim(6.0,5.0E3)
ax41.set_xlim(0.1,23.0)
ax41.set_ylim(6.0,5.0E3)

#ax.legend()
ax11.grid(True)
ax12.grid(True)
ax21.grid(True)
ax22.grid(True)
ax31.grid(True)
ax32.grid(True)
ax41.grid(True)
ax42.grid(True)

  
axes1 = [ax11,ax31,ax21,ax41]
axes2 = [ax12,ax32,ax22,ax42]

c = ['b','r','g']
lt = ['k-.','k','k--']

mass_range_v = ['12.50-12.97','12.97-13.33','13.33-13.74','13.74-15.00']
mag = 17
kind = 'CG'
angulo = '00'
bc = 1.00
ab = 1.00
al_b = 1.00
al_c = 1.00
sbc = "%.2f" % bc
sab = "%.2f" % ab
sal_b = "%.2f" % al_b
sal_c = "%.2f" % al_c

r0_v = [0.1,0.3,0.7,0.9]

for mass_range, ax, ax1, r0 in zip(mass_range_v,axes1,axes2,r0_v):

  ax.set_title(mass_range)

###  nlines = 50
###  ndir = 3
###  ym1h = np.zeros((ndir, nlines))
###  er1h = np.zeros((ndir, nlines))
###  ym2h = np.zeros((ndir, nlines))
###  er2h = np.zeros((ndir, nlines))
###  ym   = np.zeros((ndir, nlines))
###  er   = np.zeros((ndir, nlines))
###
###  suffix = mass_range+'_'+angulo+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+'%2d'%mag+'.CG.dat'
###  filename1h = '../outputs/chi/2ndpass/funcorr_1h_'+suffix
###
###  suffix = mass_range+'_'+angulo+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+'%2d'%mag+'.CG.dat'
###  filename2h = '../outputs/chi/2ndpass/funcorr_2h_'+suffix
###
###  if os.path.isfile(filename1h):
###    data = np.genfromtxt(filename1h,unpack=True)
###    x = data[0]
###    ym1h[0,:] = data[1]
###    er1h[0,:] = data[2]
###    ym1h[1,:] = data[3]
###    er1h[1,:] = data[4]
###    ym1h[2,:] = data[5]
###    er1h[2,:] = data[6]
###
###  if os.path.isfile(filename2h):
###    data = np.genfromtxt(filename2h,unpack=True)
###    x = data[0]
###    ym2h[0,:] = data[1]
###    er2h[0,:] = data[2]
###    ym2h[1,:] = data[3]
###    er2h[1,:] = data[4]
###    ym2h[2,:] = data[5]
###    er2h[2,:] = data[6]
###
###  #r0 = 0.9
###  ym2h = (erf((np.log10(x) - np.log10(r0))*3.0)*0.5 + 0.5)*ym2h
###
###  ym = ym1h + ym2h
###  #er = np.sqrt(er1h*er1h + er2h*er2h)
###
###  y = ym[0,:]
###  error = er[0,:]
###  color = 'red'
###  ax.plot(x,y, '-', lw=2,color=color)
###  ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)
###
###  y = ym[1,:]
###  error = er[1,:]
###  color = 'blue'
###  ax.plot(x,y, '-', lw=2,color=color)
###  ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)
###
###  y = ym[2,:]
###  error = er[2,:]
###  color = 'black'
###  ax.plot(x,y,'-',lw=2,color=color,label="Modelo")
###  ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)
###  
###  y = ym1h[0,:]
###  ax.plot(x,y, ':', color='black')
###  y = ym1h[1,:]
###  ax.plot(x,y, ':', color='black')
###  y = ym1h[2,:]
###  ax.plot(x,y, ':', color='black')
###  
###  y = ym2h[0,:]
###  ax.plot(x,y, ':', color='black')
###  y = ym2h[1,:]
###  ax.plot(x,y, ':', color='black')
###  y = ym2h[2,:]
###  ax.plot(x,y, ':', color='black')


#color = 'red'
#ax1.plot(x,ym[0,:]/ym[2,:], '-', lw=2,color=color)
#color = 'blue'
#ax1.plot(x,ym[1,:]/ym[2,:], '-', lw=2,color=color)
#color = 'black'
#ax1.plot(x,ym[2,:]/ym[2,:],'-',lw=2,color=color)

  modelo(mass_range,ax)

  filename = 'facu/mockscorpio/xiproy.jack_-'+'%2d'%mag+'_'+mass_range+'_20bin_100mp.dat'
  r = 0.0
  xi = xipar = xiper = 0.0
  if os.path.isfile(filename):
      r,xi,xipar,xiper=np.genfromtxt(filename,unpack=True,usecols=(0,1,3,5))
      ax.plot(r,xi,'--k',lw=1,label=r"Proyectada desde $\xi(\sigma,\pi)$")
      ax.plot(r,xipar,'--ob',lw=1,label=r"Proyectada $\perp$")
      ax.plot(r,xiper,'--or',lw=1,label=r"Proyectada $\parallel$")
  
      #cocientes
      ax1.plot(r,xi/xi,'--k')
      ax1.plot(r,xipar/xi,'--ob',lw=1)
      ax1.plot(r,xiper/xi,'--or',lw=1)

  rcParams['legend.frameon']=False

plt.show()
