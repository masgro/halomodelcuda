import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf

def lee_xi_iso(filename):

          x = np.genfromtxt(filename) ##jacknife
          xiso = x[:,0]/1000.0
          yiso = x[:,5]/x[:,6] - 1.0

          return(xiso,yiso)


def plot_xi(filename,yiso):
          c = ['b','r','g']

          x = np.genfromtxt(filename) ##jacknife
          x1 = x[:,0]/1000.0
          y1 = x[:,1]
          y2 = x[:,3]
          y3 = x[:,5]

          ax1.plot(x1,y1, 'o', color=c[0], marker='o')
          ax1.plot(x1,y2, 'o', color=c[1], marker='o')
          ax1.plot(x1,y3, 'o', color=c[2], marker='o')

          ax2.plot(x1,y1/yiso, 'o', color=c[0], marker='o')
          ax2.plot(x1,y2/yiso, 'o', color=c[1], marker='o')
          ax2.plot(x1,y3/yiso, 'o', color=c[2], marker='o')
          ax2.plot(x1,y3*0.0+1.0, '--', color='black')

def plot_hm(mmin,mmax,bc,ab,al_b,al_c):
          path = './'

          name = '45'

          sbc = "%.2f" % bc
          sab = "%.2f" % ab
          sal_b = "%.2f" % al_b
          sal_c = "%.2f" % al_c

          file = path+'funcorr_2h_'+mmin+'-'+mmax+'_'+name+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'.dat'
          x = np.genfromtxt(file)
          x1 = x[:,0]
          y1 = x[:,1]
          y2 = x[:,2]
          y3 = x[:,3]

          ax1.plot(x1, y1, '-', color='blue', linewidth = 1)
          ax1.plot(x1, y2, '-', color='red',  linewidth = 1)
          ax1.plot(x1, y3, '-', color='green',linewidth = 1)

rc('text', usetex=False)

f = plt.figure(figsize=(9,9))

gs = gridspec.GridSpec(3, 1)
gs.update(left=0.1, right=0.96, hspace=0.00)
ax1 = plt.subplot(gs[:-1,0])
ax2 = plt.subplot(gs[2,0],sharex=ax1)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xscale('log')

i = 0

m1 = ['12.00','12.50','13.00','13.50','14.00','14.50']
m2 = ['12.50','13.00','13.50','14.00','14.50','15.00']

x1hmin = [2.0,0.3,0.5,0.7,1.0,1.5]
x1hmax = [8.0,0.6,1.0,1.5,2.0,3.0]

xbest = [0.85, 0.00, 0.00, 0.00, 0.00, 0.57] 
ybest = [0.45, 0.00, 0.00, 0.00, 0.00, 0.33]

m1 = m1[i]
m2 = m2[i]
x1hmin = x1hmin[i]
x1hmax = x1hmax[i]
xbest = xbest[i]
ybest = ybest[i]

xiso, yiso = lee_xi_iso('saraswati/fccp_'+m1+'_'+m2+'_45.dat')
plot_xi('saraswati/fccp_saraswati_'+m1+'_'+m2+'_45.jacknife',yiso)

nf = 10
al_bmin = 0.1
al_bmax = 0.8
al_cmin = 0.3
al_cmax = 0.6
dal_b = (al_bmax - al_bmin)/(nf - 1)
dal_c = (al_cmax - al_cmin)/(nf - 1)

##i = 0
##while i < nf:
##  j = 0
##  while j < nf:
##    al_b = al_bmin + dal_b*i
##    al_c = al_cmin + dal_c*j
##    plot_hm(m1,m2,0.99,0.99,al_b,al_c)
##    j += 1
##  i += 1

plot_hm(m1,m2,0.99,0.99,xbest,ybest)

ax2.axvline(x= x1hmin, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax2.axvline(x= x1hmax, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax1.axvline(x= x1hmin, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)
ax1.axvline(x= x1hmax, ymin=0.0, ymax=1.0,linestyle="--",linewidth=2)


ax1.set_xlim(0.04,60.00)
ax1.set_ylim(5.0E-2,1.0E5)
ax1.set_ylabel(r'$\xi (r)$')
ax1.set_xlabel(r'$r [h^{-1} Mpc]$')

ax2.set_xlim(0.04,60.00)
ax2.set_ylim(0,2)
ax2.set_ylabel(r'$quotients$')
ax2.set_xlabel(r'$r [h^{-1} Mpc]$')

plt.show()


