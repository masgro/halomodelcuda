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

def plot_func(i,j,color):
      ii = "%1d" % i
      jj = "%02d" % j
      file = '../funcorr_'+ii+'_2h.'+jj
      data = np.genfromtxt(file)
      x = data[:,0]
      y = data[:,1]
      e = data[:,2]
      
      ax.plot(x,y, '-', color=color)
      ax.fill_between(x, y-e, y+e,alpha=0.5, edgecolor=color, facecolor=color,
                       linewidth=0, linestyle='solid', antialiased=True)


def plotcorr(kind):

      file = '../funcorr_1h_14.00-14.50_45_0.50_0.78_0.48_0.52.'+kind+'.dat'
      data = np.genfromtxt(file)

      nlines = len(data[:,0])
      y1h = np.zeros((3,nlines))
      e1h = np.zeros((3,nlines))
      y2h = np.zeros((3,nlines))
      e2h = np.zeros((3,nlines))
      y = np.zeros((3,nlines))
      e = np.zeros((3,nlines))

      x1 = data[:,0]
      y1h[0,:] = data[:,1]
      e1h[0,:] = data[:,2]
      y1h[1,:] = data[:,3]
      e1h[1,:] = data[:,4]
      y1h[2,:] = data[:,5]
      e1h[2,:] = data[:,6]

      file = '../funcorr_2h_14.00-14.50_45_0.50_0.78_0.48_0.52.'+kind+'.dat'
      data = np.genfromtxt(file)

      x2 = data[:,0]
      y2h[0,:] = data[:,1]
      e2h[0,:] = data[:,2]
      y2h[1,:] = data[:,3]
      e2h[1,:] = data[:,4]
      y2h[2,:] = data[:,5]
      e2h[2,:] = data[:,6]

      e = np.power(e1h,2) + np.power(e2h,2)
      e = np.sqrt(e)

      x = x1

      y = y1h[0,:] + y2h[0,:]
      error = e[0,:]
      color = 'red'
      ax.plot(x,y,linewidth=1.5, linestyle="-", color=color)
      ax.fill_between(x, y-error, y+error,
                      alpha=0.5, edgecolor=color, facecolor=color,
                      linewidth=0, linestyle='solid', antialiased=True)

      y = y1h[1,:] + y2h[1,:]
      error = e[1,:]
      color = 'blue'
      ax.plot(x,y,linewidth=1.5, linestyle="-", color=color)
      ax.fill_between(x, y-error, y+error,
                      alpha=0.5, edgecolor=color, facecolor=color,
                      linewidth=0, linestyle='solid', antialiased=True)

      y = y1h[2,:] + y2h[2,:]
      error = e[2,:]
      color = 'green'
      ax.plot(x,y, linewidth=1.5, linestyle="-", color=color)
      ax.fill_between(x, y-error, y+error,
                      alpha=0.5, edgecolor=color, facecolor=color,
                      linewidth=0, linestyle='solid', antialiased=True)


rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

f = plt.figure(figsize=(8,8))

ax = f.add_subplot(111)

ax.set_xscale('log')
ax.set_yscale('log')

plotcorr('HM')
plotcorr('CG')

ax.set_xlim(0.1,50.0)
ax.set_ylim(1.0E-1,1.0E4)
ax.set_ylabel(r'$\xi (r)$')
ax.set_xlabel(r'$r [h^{-1} Mpc]$')

ax.grid(True)

plt.show()
