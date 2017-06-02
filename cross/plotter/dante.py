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
          ax1.plot(x1,y4, 'o', color=c[2], marker='o')

          ax2.plot(x1,y1/y4, 'o', color=c[0], marker='o')
          ax2.plot(x1,y2/y4, 'o', color=c[1], marker='o')
          ax2.plot(x1,y3/y4, 'o', color=c[2], marker='o')

def plot_hm_dante(mmin,mmax,bc,ab,al_b,al_c):
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
          y1 = x[:,1]
          y2 = x[:,2]
          y3 = x[:,3]

          #file = 'funcorr_2h_'+mmin+'-'+mmax+'_'+name+'_'+sbc+'_'+sab+'_1.000_0.000.dat'
          #x = np.genfromtxt(file)
          #y1 += al_b*x[:,1]
          #y2 += al_b*x[:,2]
          #y3 += al_b*x[:,3]

          #file = 'funcorr_2h_'+mmin+'-'+mmax+'_'+name+'_'+sbc+'_'+sab+'_0.000_3.000.dat'
          #x = np.genfromtxt(file)
          #y1 += al_c*x[:,1]
          #y2 += al_c*x[:,2]
          #y3 += al_c*x[:,3]

          #norma = 1.0 + al_b*1.2568 + 0.5*1.2568*al_c
          #norma = 1.0 + al_b*1.2568 + 1.2568*al_c
          norma = 1.0 + al_b*1.2568 + (log(4.0) - 1.0)*1.2568*al_c
          norma = 1./norma

          y2h[:,0] = y1#*norma*0.81368201;
          #y2h[:,1] = y2*norma*0.81368201;
          #y2h[:,2] = y3*norma*0.81368201;
          
          lt = ['k-.','k','k--']
          
          y = y2h[:,0]
          ax1.plot(x1, y, lt[0], color='red',linewidth = 2)

          #y = y2h[:,1]
          #ax1.plot(x1, y, lt[1], color='green',linewidth = 2)
  
          #y = y2h[:,2]
          #ax1.plot(x1, y, lt[2], color='blue',linewidth = 2)

def plot_hm(mmin,mmax,bc,ab,al_b,al_c):
          c = ['b','r','g']
          
          nlines = 100
          name = 'cono'

          y2h = np.zeros((nlines,3))
          
          sbc = "%.2f" % bc
          sab = "%.2f" % ab
          sal_b = "%.2f" % al_b
          sal_c = "%.2f" % al_c

          file = 'funcorr_1h_'+mmin+'-'+mmax+'_'+name+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'.dat'
          x = np.genfromtxt(file)

          x1 = x[:,0]
          y1 = x[:,1]
          y2 = x[:,2]
          y3 = x[:,3]

          y2h[:,0] = y1
          y2h[:,1] = y2
          y2h[:,2] = y3
          
          lt = ['-','-','-']
          
          y = y2h[:,0]
          ax1.plot(x1, y, lt[0], color='red',linewidth = 1)

          y = y2h[:,1]
          ax1.plot(x1, y, lt[1], color='green',linewidth = 1)
  
          y = y2h[:,2]
          ax1.plot(x1, y, lt[2], color='blue',linewidth = 1)
  

rc('text', usetex=False)

f = plt.figure(figsize=(9,9))

gs = gridspec.GridSpec(3, 1)
gs.update(left=0.1, right=0.96, hspace=0.00)
ax1 = plt.subplot(gs[:-1,0])
ax2 = plt.subplot(gs[2,0],sharex=ax1)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xscale('log')

m1 = '14.50'
m2 = '15.00'

plot_xi('saraswati/fccp_12.00_12.50_45.dat')

plot_xi('saraswati/fccp_14.50_15.00_45.dat')

#12.00 - 12.50
#x1hmin = 0.2
#x1hmax = 0.4
#x2hmin = 2.0
#x2hmax = 8.0

#12.50 - 13.00
#x1hmin = 0.3
#x1hmax = 0.6
#x2hmin = 2.0
#x2hmax = 9.0

#13.00 - 13.50
#x1hmin = 0.5
#x1hmax = 1.0
#x2hmin = 3.0
#x2hmax = 10.0

#13.50 - 14.00
#x1hmin = 0.7
#x1hmax = 1.5
#x2hmin = 4.0
#x2hmax = 10.0

#14.00 - 14.50
#x1hmin = 1.0
#x1hmax = 2.0
#x2hmin = 6.0
#x2hmax = 15.0

#14.50 - 15.00
x1hmin = 1.5
x1hmax = 3.0
x2hmin = 7.0
x2hmax = 15.0

ax2.axvline(x= x1hmin, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax2.axvline(x= x1hmax, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax2.axvline(x= x2hmin, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax2.axvline(x= x2hmax, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax1.axvline(x= x1hmin, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax1.axvline(x= x1hmax, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax1.axvline(x= x2hmin, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax1.axvline(x= x2hmax, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)

plot_hm_dante('12.00','12.50',0.79,0.89,0.0,0.0)

plot_hm_dante('14.50','15.00',0.53,0.68,0.0,0.0)
#plot_hm_dante(m1,m2,0.7,0.8,14.24,5.384) #con termino de acople
#plot_hm_dante(m1,m2,0.7,0.8,0.313,0.202) #desacopladas

#r = 2.182
#alpha = 0.323
#plot_hm_dante(m1,m2,0.7,0.8,r*cos(alpha),r*sin(alpha)) #con termino de acople doble
#

nf = 20

#al_bmin = 0.0
#al_bmax = 1.0
#al_cmin = 0.0
#al_cmax = 1.0
#dal_b = (al_bmax - al_bmin)/(nf - 1)
#dal_c = (al_cmax - al_cmin)/(nf - 1)

#bcmin = 0.0
#bcmax = 1.0
#abmin = 0.0
#abmax = 1.0
#dbc = (bcmax - bcmin)/(nf - 1)
#dab = (abmax - abmin)/(nf - 1)
#
#i = 0
#while i < nf:
#  j = 0
#  while j < nf:
#        #al_b = al_bmin + dal_b*i
#        #al_c = al_cmin + dal_c*j
#        #al_b = np.exp(al_b)
#        #al_c = al_b/2.525
#        #plot_hm_dante(m1,m2,0.7,0.8,np.abs(al_b),np.abs(al_c))
#
#        bc = bcmin + dbc*i
#        ab = abmin + dab*j
#        plot_hm(m1,m2,bc,ab,0.0,0.0)
#        
#        j += 1
#  i += 1

ax1.set_xlim(0.04,60.00)
ax1.set_ylim(5.0E-2,1.0E5)
ax1.set_ylabel(r'$\xi (r)$')
ax1.set_xlabel(r'$r [h^{-1} Mpc]$')

ax2.set_xlim(0.04,60.00)
ax2.set_ylim(0,2)
ax1.set_ylabel(r'$\xi (r)$')
ax1.set_xlabel(r'$r [h^{-1} Mpc]$')

plt.show()
