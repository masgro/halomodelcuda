import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf

def plot_xi(filename):
          c = ['b','r','g']
          ##x = np.genfromtxt(filename) ##jacknife
          ##x1 = x[:,0]/1000.0
          ##y1 = x[:,1]
          ##y2 = x[:,3]
          ##y3 = x[:,5]

          x = np.genfromtxt(filename) ##jacknife
          x1 = x[:,0]/1000.0
          y1 = x[:,1]/x[:,4] - 1.0
          y2 = x[:,2]/x[:,4] - 1.0
          y3 = x[:,3]/x[:,4] - 1.0
          y4 = x[:,5]/x[:,6] - 1.0

          #ax1.plot(x1,y1, 'o', color=c[0], marker='o')
          #ax1.plot(x1,y2, 'o', color=c[1], marker='o')
          #ax1.plot(x1,y3, 'o', color=c[2], marker='o')
          ax1.plot(x1,y4, 'o', color='black', marker='o')

def plot_hm_dante(mmin,mmax,bc,ab,al_b,al_c,norm):
          c = ['b','r','g']
          
          nlines = 50
          name = '90'

          y2h = np.zeros((nlines,3))
          
          sbc = "%.2f" % bc
          sab = "%.2f" % ab
          sal_b = "%.2f" % al_b
          sal_c = "%.2f" % al_c

          file = 'funcorr_2h_'+mmin+'-'+mmax+'_'+name+'_'+sbc+'_'+sab+'_0.00_0.00.dat'
          x = np.genfromtxt(file)
          x1 = x[:,0]
          y1 = x[:,1]*norm
          y2 = x[:,2]*norm
          y3 = x[:,3]*norm

          ax1.plot(x1, y1, '-', color='red',linewidth = 2)

rc('text', usetex=False)

f = plt.figure(figsize=(9,9))

ax1 = plt.subplot(111)

ax1.set_xscale('log')
ax1.set_yscale('log')

m1 = ['12.00','12.50','13.00','13.50','14.00','14.50']
m2 = ['12.50','13.00','13.50','14.00','14.50','15.00']

x1hmin = [0.2,0.3,0.5,0.7,1.0,1.5]
x1hmax = [0.4,0.6,1.0,1.5,2.0,3.0]

x2hmin = [2.0,2.0,3.0,4.0,6.0,7.0]
x2hmax = [8.0,9.0,10.0,10.0,15.0,15.0]

xbest = [0.79, 0.79, 0.68, 0.58, 0.58, 0.53] 
ybest = [0.89, 0.84, 0.89, 0.84, 0.74, 0.68]


norm = [0.778992,0.738978,0.832521,0.929103,1.006819,1.172130]

i = 0
while i < 6:
  mmin = m1[i]
  mmax = m2[i]
  bc = xbest[i]
  ab = ybest[i]
  plot_xi('saraswati/fccp_'+mmin+'_'+mmax+'_45.dat')
  plot_hm_dante(mmin,mmax,bc,ab,0.0,0.0,norm[i])
  i = i + 1


ax1.set_xlim(0.04,60.00)
ax1.set_ylim(5.0E-2,1.0E5)
ax1.set_ylabel(r'$\xi (r)$')
ax1.set_xlabel(r'$r [h^{-1} Mpc]$')

plt.show()
