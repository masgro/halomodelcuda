import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf

def xilin(x):
          "Funcion de correlacion lineal"
          result  = 0.804728
          result -= 0.963419*x
          result -= 0.408048*np.power(x,2)
          result += 0.056507*np.power(x,3)
          result += 0.096629*np.power(x,4)
          result -= 0.049933*np.power(x,5)
          result -= 0.032321*np.power(x,6)
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
mmin = '12.00'
mmax = '12.50'

name = 'cono'
angulo = '45'
nfiles = 5
nlines = 50

ym = np.zeros((3, nfiles, nlines))
y1h = np.zeros((3, nfiles, nlines))
y2h = np.zeros((3, nfiles, nlines))

factor = 0.9

r0 = 0.5

b = 0
while b < nfiles:
  a = 0
  while a < 3:
    #file = 'data/funcorr_'+str(a)+'_1h_'+mmin+'-'+mmax+'_'+name+'_'+angulo+'_99.dat'
    #x = np.genfromtxt(file)
    #x1 = x[:,0]
    #y1 = x[:,1]*factor
  
    #file = 'data/funcorr_'+str(a)+'_2h_'+mmin+'-'+mmax+'_'+name+'_'+angulo+'_'+str(b).zfill(2)+'.dat'
    file = 'data/funcorr_'+str(a)+'_2h.'+str(b).zfill(2)
    x = np.genfromtxt(file)
    x2 = x[:,0]
    y2 = x[:,1]*factor*factor*(erf((np.log10(x2) - np.log10(r0))*3.0)*0.5 + 0.5)*0.7

    x1 = x2
    #ym[a,b,:] = (y1 + y2)
    ym[a,b,:] = y2
    #y1h[a,b,:] = y1
    y2h[a,b,:] = y2

    a += 1
  b += 1

lt = ['k-.','k','k--']

yy = np.zeros((3, nlines))

a = 0
while a < 3:
  y = ym[a,:,:]
  y.reshape(nfiles,nlines)
  y = np.median(y, axis=0)
  yy[a,:] = y
  ax1.plot(x1,y, lt[a], color='black',linewidth = 2)

  y = y1h[a,:,:]
  y.reshape(nfiles,nlines)
  y = np.median(y, axis=0)
  ax1.plot(x1,y, ':', color='black',linewidth = 2)

  y = y2h[a,:,:]
  y.reshape(nfiles,nlines)
  y = np.median(y, axis=0)
  ax1.plot(x1,y, ':', color='black',linewidth = 2)

  a += 1

file = 'saraswati/fccp_'+mmin+'_'+mmax+'_'+angulo+'.dat'
x = np.genfromtxt(file)
x1 = x[:,0]/1000.0
y1 = xi(x[:,1],x[:,4])
y2 = xi(x[:,2],x[:,4])
y3 = xi(x[:,3],x[:,4])
y4 = xi(x[:,5],x[:,6])

ax1.plot(x1,y1, 'o', color='black', marker='o', mfc='None')
ax1.plot(x1,y2, 'o', color='black', marker='d', mfc='None')
ax1.plot(x1,y3, 'o', color='black', marker='s', mfc='None')
ax1.plot(x1,y4, 'o', color='black', marker='o')

ax2.plot(x1,y1/y4, 'o', color='black', marker='o', mfc='None')
ax2.plot(x1,y2/y4, 'o', color='black', marker='d', mfc='None')
ax2.plot(x1,y3/y4, 'o', color='black', marker='s', mfc='None')

name = 'cono'
nfiles = 5
nlines = 100
angulo = '90'
ym = np.zeros((3, nfiles, nlines))
y1h = np.zeros((3, nfiles, nlines))
y2h = np.zeros((3, nfiles, nlines))

b = 0
while b < nfiles:
  a = 0
  while a < 3:
    file = 'data/funcorr_'+str(a)+'_1h_'+mmin+'-'+mmax+'_'+name+'_'+angulo+'_'+str(b).zfill(2)+'.dat'
    x = np.genfromtxt(file)
    x1 = x[:,0]
    y1 = x[:,1]*factor
  
    file = 'data/funcorr_'+str(a)+'_2h_'+mmin+'-'+mmax+'_'+name+'_'+angulo+'_'+str(b).zfill(2)+'.dat'
    x = np.genfromtxt(file)
    x2 = x[:,0]
    y2 = x[:,1]*factor*factor*(erf((np.log10(x2) - np.log10(r0))*3.0)*0.5 + 0.5)*0.7

    ym[a,b,:] = (y1 + y2)
    y1h[a,b,:] = y1
    y2h[a,b,:] = y2

    a += 1
  b += 1

y = ym[0,:,:]
y.reshape(nfiles,nlines)
y = np.median(y, axis=0)
ax1.plot(x1,y,'k', color='blue')

y = ym[1,:,:]
y.reshape(nfiles,nlines)
y = np.median(y, axis=0)
ax1.plot(x1,y,'k', color='blue')

y = ym[2,:,:]
y.reshape(nfiles,nlines)
y = np.median(y, axis=0)
ax1.plot(x1,y,'k', color='blue')

lt = ['k-.','k','k--']
ax2.plot(x1,yy[0,:]/y, lt[0], color='black')
ax2.plot(x1,yy[1,:]/y, lt[1], color='black')
ax2.plot(x1,yy[2,:]/y, lt[2], color='black')

ax1.set_xlim(0.1,40.0)
ax1.set_ylim(5.0E-2,1.0E4)
ax1.set_ylabel(r'$\xi (r)$')

#linea horizotal en 1
ax2.plot(x1,x1*0.0 + 1.0, ':', color='black')

ax2.set_xlim(0.1,40.0)
ax2.set_ylim(0.0,2.0)
ax2.set_xlabel(r'$r [Mpc \  h^{-1}]$')
ax2.set_ylabel(r'$quotients$')

plt.show()

#fileout = 'cono.ps'
#plt.savefig(fileout)
