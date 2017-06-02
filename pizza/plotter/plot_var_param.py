import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf
from scipy import interpolate
import os.path
from scipy.interpolate import interp1d
from scipy.optimize import minimize

ndir = 3
nlines = 50

angulo = '45'
mmin = 13.33
mmax = 13.74
MAG = 17


def myfunc(par,r,xi,xipar,xiper,x,y1h,y2h):

  f1h0 = interp1d(x,y1h[0,:])
  f1h1 = interp1d(x,y1h[1,:])
  f1h2 = interp1d(x,y1h[2,:])
  f2h0 = interp1d(x,y2h[0,:])
  f2h1 = interp1d(x,y2h[1,:])
  f2h2 = interp1d(x,y2h[2,:])

  s0 = s1 = s2 = 0.0

  rmin = 0.4
  rmax = 11.0
  indx = np.where((r > rmin) & (r < rmax))
  xx = r[indx]

  yy = xi[indx]
  ytmp2h = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*f2h2(xx)
  ytmp1h = (erf((np.log10(par[0]) - np.log10(xx))*par[1])*0.5 + 0.5)*f1h2(xx)
  ytmp = ytmp1h + ytmp2h 
  s2 = np.sum(np.power(yy-ytmp,2.))

  #yy = xipar[indx]
  #ytmp2h = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*f2h1(xx)
  #ytmp1h = (erf((np.log10(par[0]) - np.log10(xx))*par[1])*0.5 + 0.5)*f1h1(xx)
  #ytmp = ytmp1h + ytmp2h 
  #s1 = np.sum(np.power(yy-ytmp,2.))

  #yy = xiper[indx]
  #ytmp2h = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*f2h0(xx)
  #ytmp1h = (erf((np.log10(par[0]) - np.log10(xx))*par[1])*0.5 + 0.5)*f1h0(xx)
  #ytmp = ytmp1h + ytmp2h 
  #s0 = np.sum(np.power(yy-ytmp,2.))

  return(s0+s1+s2)

rc('text', usetex=False)
rcParams.update({'font.size': 8})

f = plt.figure(figsize=(8,8))

gs = gridspec.GridSpec(1,1)
gs.update(left=0.1, right=0.96, hspace=0.00, wspace=0.00)
ax = plt.subplot(gs[0])

ax.set_xscale('log')
ax.set_yscale('log')

filename = 'facu/mockscorpio/xiproy.jack_-'+'%2d'%MAG+'_'+'%.2f'%mmin+'-'+'%.2f'%mmax+'_20bin_100mp.dat'
if os.path.isfile(filename):
  r,xi,xipar,xiper=np.genfromtxt(filename,unpack=True,usecols=(0,1,3,5))
  ax.plot(r,xi,'--ok',lw=1,label=r"Proyectada desde $\xi(\sigma,\pi)$")
  ax.plot(r,xipar,'--ob',lw=1,label=r"Proyectada $\perp$")
  ax.plot(r,xiper,'--or',lw=1,label=r"Proyectada $\parallel$")

    #cocientes
#    ax1.plot(r,xi/xi,'--ok')
#    ax1.plot(r,xipar/xi,'--ob',lw=1)
#    ax1.plot(r,xiper/xi,'--or',lw=1)


def plot_hm_tipo(mmin,mmax,angle,bc,ab,al_b,al_c,mag,tipo):

  y1h = np.zeros((ndir, nlines))
  y2h = np.zeros((ndir, nlines))

  smmin = "%.02f" % mmin
  smmax = "%.02f" % mmax
  sangle = "%02d" % angle
  sbc = "%.2f" % bc
  sab = "%.2f" % ab
  sal_b = "%.2f" % al_b #1.80
  sal_c = "%.2f" % al_c #1.80
  smag = "%2d" % mag
  file = '../outputs/funcorr_1h_'+smmin+'-'+smmax+'_'+sangle+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+smag+'.CG.dat'

  if os.path.isfile(file):
    data = np.genfromtxt(file)
    x = data[:,0]
    y1h[0] = data[:,1]#*1.2#*1.9#
    y1h[1] = data[:,3]#*1.2#*1.9#
    y1h[2] = data[:,5]#*1.2#*1.9#
  else :
    print file
    y1h = 0.0

  sbc = "%.2f" % 1.00
  sab = "%.2f" % 1.00
  sal_b = "%.2f" % al_b
  sal_c = "%.2f" % al_c

  file = '../outputs/funcorr_2h_'+smmin+'-'+smmax+'_'+sangle+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+smag+'.CG.dat'
  if os.path.isfile(file):
    data = np.genfromtxt(file)
    x = data[:,0]
    y2h[0] = data[:,1]#*1.15
    y2h[1] = data[:,3]#*1.15
    y2h[2] = data[:,5]#*1.15
  else :
    print file
    y2h = 0.0


  #x0 = [1.6, 7.7]
  #res = minimize(myfunc, x0, args=(r,xi,xipar,xiper,x,y1h,y2h),  method='Nelder-Mead',bounds=[[0.3,23.0],[1.0,7.0]])
  #x0 = res.x


  #y2h *= (erf((np.log10(x) - np.log10(x0[0]))*x0[1])*0.5+0.5)
  #y1h *= (erf((np.log10(x0[0]) - np.log10(x))*x0[1])*0.5+0.5)

  return x,y1h+y2h

Avec = np.linspace(0.1,0.9,2)
Bvec = np.linspace(0.1,0.9,2)

Cvec = np.linspace(0.1,1.0,9)
Dvec = np.linspace(0.1,1.0,9)

Colorvec = ['red','blue']

for A,color in zip(Avec,Colorvec):
    for B in Bvec:
        #for C in Cvec :
            #for D in Dvec :

                #B = 0.1
                C = 1.8
                D = 1.8
                x,y = plot_hm_tipo(mmin,mmax,45,A,B,C,D,MAG,1)
                #ax.plot(x,y[0],'-',color='red',linewidth = 1)
                #ax.plot(x,y[1],'-',color='blue',linewidth = 1)
                ax.plot(x,y[2],'-',color=color,linewidth = 1)

Colorvec = ['green','black']
for A,color in zip(Avec,Colorvec):
    for B in Bvec:
        #for C in Cvec :
            #for D in Dvec :

                #B = 0.1
                C = 1.9
                D = 1.9
                x,y = plot_hm_tipo(mmin,mmax,45,A,B,C,D,MAG,1)
                #ax.plot(x,y[0],'-',color='red',linewidth = 1)
                #ax.plot(x,y[1],'-',color='blue',linewidth = 1)
                ax.plot(x,y[2],'-',color=color,linewidth = 1)


filename = 'facu/mockscorpio/xiproy.jack_-'+'%2d'%MAG+'_'+'%.2f'%mmin+'-'+'%.2f'%mmax+'_20bin_100mp.dat'
if os.path.isfile(filename):
  r,xi,xipar,xiper=np.genfromtxt(filename,unpack=True,usecols=(0,1,3,5))
  ax.plot(r,xi,'--ok',lw=1,label=r"Proyectada desde $\xi(\sigma,\pi)$")
  ax.plot(r,xipar,'--ob',lw=1,label=r"Proyectada $\perp$")
  ax.plot(r,xiper,'--or',lw=1,label=r"Proyectada $\parallel$")

    #cocientes
#    ax1.plot(r,xi/xi,'--ok')
#    ax1.plot(r,xipar/xi,'--ob',lw=1)
#    ax1.plot(r,xiper/xi,'--or',lw=1)


majorFormatter = FormatStrFormatter('%.1f')

ax.set_xlim(0.1,30.00)
ax.set_ylim(6.0,5.0E3)
ax.set_ylabel(r'$\xi (r)$')


plt.show()
