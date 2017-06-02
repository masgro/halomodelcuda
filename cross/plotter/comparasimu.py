import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec


def xilin(x):
          "Funcion de correlacion lineal"
          #result = 6.951502249814009e-01
          #result = result - 8.559727903289444e-01*x
          #result = result - 2.263312579328268e-01*np.power(x,2)
          #result = result - 4.379787572393435e-02*np.power(x,3)
          #result = result - 9.014241771627268e-03*np.power(x,4)
          #result = result - 1.436052315861557e-02*np.power(x,5)
          #result = result - 5.332102016644595e-03*np.power(x,6)

          result =  0.6806686591027654
          result = result - 0.8151080596170317*x
          result = result - 0.1694839355035661*np.power(x,2)
          result = result - 0.0759901334068769*np.power(x,3)
          result = result - 0.0535756512914830*np.power(x,4)


          return result

def xi(x,y):
          "Computa la funcion de correlacion"
          result = x/y - 1.0
          return result

rc('text', usetex=False)

f = plt.figure(figsize=(9,12))

gs = gridspec.GridSpec(3, 1)
gs.update(left=0.1, right=0.96, hspace=0.0)

ax1 = plt.subplot(gs[:-1,0])
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2 = plt.subplot(gs[2,0],sharex=ax1)
ax2.set_xscale('log')

c = ['b','r','g']
mmin = '12.0'
mmax = '12.5'

angulo = '45'

file = 'fccp_saraswati_12.00_12.50_'+angulo+'.dat'
x = np.genfromtxt(file)
x1 = x[:,0]/1000.0
y1 = xi(x[:,1],x[:,4])
y2 = xi(x[:,2],x[:,4])
y3 = xi(x[:,3],x[:,4])
y4 = xi(x[:,5],x[:,6])

#ax1.plot(x1,y1, 'k', color='black')
#ax1.plot(x1,y2, 'k', color='black')
#ax1.plot(x1,y3, 'k', color='black')
ax1.plot(x1,y4, 'o', color='black', marker='d')

ax2.plot(x1,y1/y4, 'k', color='black')
ax2.plot(x1,y2/y4, 'k', color='black')
ax2.plot(x1,y3/y4, 'k', color='black')

ax1.set_xlim(0.1,40.0)
ax1.set_ylim(5.0E-2,1.0E4)
ax1.set_ylabel(r'$\xi (r)$')
ax1.get_xaxis().set_ticklabels([])

ax2.set_xlim(0.1,40.0)
ax2.set_ylim(0.0,2.0)
ax2.set_xlabel(r'$r [Mpc h^{-1}]$')
ax2.set_ylabel(r'$cociente$')

plt.show()

#fileout = 'cono.ps'
#plt.savefig(fileout)
