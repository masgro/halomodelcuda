#!/usr/bin/env python

import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf,gamma
import os.path
from scipy import integrate
from matplotlib.ticker import NullFormatter
from scipy.interpolate import interp1d
from scipy.optimize import minimize

#plt.style.use('grayscale')

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

f = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(4,1,)
ax = f.add_subplot(gs[0:3,0])
ax1 = f.add_subplot(gs[3,0],sharex=ax)

gs.update(hspace=0.0)
plt.setp(ax.get_xticklabels(), visible=False)

nlines = 50
ndir = 3
ym1h = np.zeros((ndir, nlines))
er1h = np.zeros((ndir, nlines))
ym2h = np.zeros((ndir, nlines))
er2h = np.zeros((ndir, nlines))
ym   = np.zeros((ndir, nlines))
er   = np.zeros((ndir, nlines))

mag = 17
#mass_range = '12.50-12.97'
#mass_range = '12.97-13.33'
mass_range = '13.33-13.74'
#mass_range = '13.74-15.00'
kind = 'CG'

angulo = '45'

bc = 1.00
ab = 1.00
al_b = 1.50
al_c = 1.50

#bc = 0.50
#ab = 0.20
#al_b = 0.55
#al_c = 0.30

sbc = "%.2f" % bc
sab = "%.2f" % ab
sal_b = "%.2f" % al_b
sal_c = "%.2f" % al_c

## PLOT MEDIDA #################################################################
#filename = 'facu/mockscorpio/xiproy_-'+'%2d'%mag+'_'+mass_range+'_20bin_100mp.dat'
#filename = 'facu/mockscorpio/xiproy_-'+'%2d'%mag+'_'+mass_range+'.dat'
#filename = 'facu/mockscorpio/xiproy_-'+'%2d'%mag+'_'+mass_range+'_20bin_100mp.mas11.dat'
filename = 'facu/mockscorpio/xiproy.jack_-'+'%2d'%mag+'_'+mass_range+'_20bin_100mp.dat'
if os.path.isfile(filename):
    r,xi,xipar,xiper=np.genfromtxt(filename,unpack=True,usecols=(0,1,3,5))
    ax.plot(r,xi,'--ok',lw=1,label=r"Proyectada desde $\xi(\sigma,\pi)$")
    ax.plot(r,xipar,'--ob',lw=1,label=r"Proyectada $\perp$")
    ax.plot(r,xiper,'--or',lw=1,label=r"Proyectada $\parallel$")

    #cocientes
    ax1.plot(r,xi/xi,'--ok')
    ax1.plot(r,xipar/xi,'--ob',lw=1)
    ax1.plot(r,xiper/xi,'--or',lw=1)
################################################################################

## PLOT MODELO #################################################################
folder = '../outputs/'
suffix = mass_range+'_'+angulo+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+'%2d'%mag+'.'+kind+'.dat'
print suffix

filename1h = folder+'funcorr_1h_'+suffix
filename2h = folder+'funcorr_2h_'+suffix

if os.path.isfile(filename1h):
    print 'reading file : ',filename1h
    data = np.genfromtxt(filename1h,unpack=True)
    x = data[0]
    ym1h[0,:] = data[1]#*1.9#*1.2
    er1h[0,:] = data[2]#*1.9#*1.2
    ym1h[1,:] = data[3]#*1.9#*1.2
    er1h[1,:] = data[4]#*1.9#*1.2
    ym1h[2,:] = data[5]#*1.9#*1.2
    er1h[2,:] = data[6]#*1.9#*1.2

if os.path.isfile(filename2h):
    print 'reading file : ',filename2h
    data = np.genfromtxt(filename2h,unpack=True)
    x = data[0]
    ym2h[0,:] = data[1]#*1.15
    er2h[0,:] = data[2]#*1.15
    ym2h[1,:] = data[3]#*1.15
    er2h[1,:] = data[4]#*1.15
    ym2h[2,:] = data[5]#*1.15
    er2h[2,:] = data[6]#*1.15

#ym1h = ym1h/9.491697e-01*1.5 #0.70 0.50
#ym1h = ym1h/8.060124e-1*1.2 #0.50 0.40
#ym1h = ym1h/1.077420e+00*1.2 #0.99 0.99

f1h0 = interp1d(x,ym1h[0,:])
f1h1 = interp1d(x,ym1h[1,:])
f1h2 = interp1d(x,ym1h[2,:])
f2h0 = interp1d(x,ym2h[0,:])
f2h1 = interp1d(x,ym2h[1,:])
f2h2 = interp1d(x,ym2h[2,:])

escala = np.zeros((3,2))
sigma = np.zeros((3,2))

if(mass_range == '12.50-12.97'):
    escala[0,1] = 0.1        ; sigma[0,1] = 0.96548314
    escala[0,0] = 0.32641641 ; sigma[0,0] = 6.99992726
    escala[1,1] = 0.55045432 ; sigma[1,1] = 4.02410762
    escala[1,0] = 0.93740779 ; sigma[1,0] = 0.74814675
    escala[2,1] = 0.23862464 ; sigma[2,1] = 1.73605431
    escala[2,0] = 0.32167312 ; sigma[2,0] = 6.58714271

if(mass_range == '12.97-13.33'):
    escala[0,1] = 0.53343258; sigma[0,1] = 3.29455895
    escala[0,0] = 0.41179681; sigma[0,0] = 3.47114219
    escala[1,1] = 0.44917009; sigma[1,1] = 3.12758295
    escala[1,0] = 0.35662938; sigma[1,0] = 4.23206749
    escala[2,1] = 0.47374043; sigma[2,1] = 2.25528555
    escala[2,0] = 0.4137295 ; sigma[2,0] = 4.0970165 

if(mass_range == '13.33-13.74'):
    #escala[0,1] = 0.50; sigma[0,1] = 1.56
    #escala[0,0] = 0.66; sigma[0,0] = 5.87
    #escala[1,1] = 1.12; sigma[1,1] = 3.46
    #escala[1,0] = 0.89; sigma[1,0] = 1.56
    #escala[2,1] = 1.07; sigma[2,1] = 4.51
    #escala[2,0] = 0.76; sigma[2,0] = 3.87

    escala[0,1] = 0.25; sigma[0,1] = 10. 
    escala[0,0] = 0.40; sigma[0,0] = 2.69
    escala[1,1] = 0.25; sigma[1,1] = 10. 
    escala[1,0] = 0.40; sigma[1,0] = 2.69
    escala[2,1] = 0.25; sigma[2,1] = 10. 
    escala[2,0] = 0.40; sigma[2,0] = 2.69

if(mass_range == '13.74-15.00'):
    escala[0,1] = 0.06709902; sigma[0,1] = 0.65937736
    escala[0,0] = 1.22901972; sigma[0,0] = 4.2920803
    escala[1,1] = 1.69025408; sigma[1,1] = 3.47046722
    escala[1,0] = 1.42679342; sigma[1,0] = 3.8874773
    escala[2,1] = 0.77660862; sigma[2,1] = 1.88762472
    escala[2,0] = 1.23729639; sigma[2,0] = 4.38157754

#for i in range(3):
#    ym1h[i] *= (erf((np.log10(escala[i,0]) - np.log10(x))*sigma[i,0])*0.5 + 0.5)
#    ym2h[i] *= (erf((np.log10(x) - np.log10(escala[i,1]))*sigma[i,1])*0.5 + 0.5)

ym = ym1h + ym2h
er = np.sqrt(er1h*er1h + er2h*er2h)

y = ym[0,:]
error = er[0,:]
color = 'red'
ax.plot(x,y, '-', lw=2,color=color)
#ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)

y = ym[1,:]
error = er[1,:]
color = 'blue'
ax.plot(x,y, '-', lw=2,color=color)
#ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)

y = ym[2,:]
error = er[2,:]
color = 'black'
ax.plot(x,y,'-',lw=2,color=color,label="Modelo")
#ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)

#y = ym1h[0,:]
#ax.plot(x,y, ':', color='black')
#y = ym1h[1,:]
#ax.plot(x,y, ':', color='black')
y = ym1h[2,:]
ax.plot(x,y, ':', color='black')

#y = ym2h[0,:]
#ax.plot(x,y, ':', color='black')
#y = ym2h[1,:]
#ax.plot(x,y, ':', color='black')
y = ym2h[2,:]
ax.plot(x,y, ':', color='black')

#cocientes
interp = interp1d(x,ym[0,:])
color = 'red'
ax1.plot(r,interp(r)/xi, '-', lw=2,color=color)
color = 'blue'
interp = interp1d(x,ym[1,:])
ax1.plot(r,interp(r)/xi, '-', lw=2,color=color)
color = 'black'
interp = interp1d(x,ym[2,:])
ax1.plot(r,interp(r)/xi, '-', lw=2,color=color)
#############################################################################################


## DETERMINA MEJORES PARAMETROS PARA EL CORTE EN 2HALOS #####################################
def myfunc(par,xt,yt,func1,func2):
    s = 0.0

    rmin = 0.11
    rmax = 6.0
    indx = np.where((xt > rmin) & (xt < rmax))
    xx = xt[indx]
    yy = yt[indx]

    yy2 = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*func2(xx)
    yy1 = (erf((np.log10(par[2]) - np.log10(xx))*par[3])*0.5 + 0.5)*func1(xx)
    yy0 = yy1 + yy2 
    s = np.sum(np.power(np.log10(yy)-np.log10(yy0),2.))

    #yy = xi[indx]
    #ytmp2h = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*f2h2(xx)
    #ytmp1h = (erf((np.log10(par[0]) - np.log10(xx))*par[1])*0.5 + 0.5)*f1h2(xx)
    #ytmp = ytmp1h + ytmp2h 
    #s2 = np.sum(np.power(np.log10(yy)-np.log10(ytmp),2.))

    #yy = xipar[indx]
    #ytmp2h = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*f2h1(xx)
    #ytmp1h = (erf((np.log10(par[0]) - np.log10(xx))*par[1])*0.5 + 0.5)*f1h1(xx)
    #ytmp = ytmp1h + ytmp2h 
    #s1 = np.sum(np.power(np.log10(yy)-np.log10(ytmp),2.))

    #yy = xiper[indx]
    #ytmp2h = (erf((np.log10(xx) - np.log10(par[0]))*par[1])*0.5 + 0.5)*f2h0(xx)
    #ytmp1h = (erf((np.log10(par[0]) - np.log10(xx))*par[1])*0.5 + 0.5)*f1h0(xx)
    #ytmp = ytmp1h + ytmp2h 
    #s0 = np.sum(np.power(np.log10(yy)-np.log10(ytmp),2.))

    return(s)

#bounds = [[0.3,3.0],[1.0,7.0]]
#x0 = [1.6, 7.7]
#res = minimize(myfunc, x0, method='Nelder-Mead',bounds=bounds)

bounds = [[0.1,10.0],[0.1,10.0],[0.1,10.0],[0.1,10.0]]
x0 = [0.5, 3.0, 0.5, 3.0]

res = minimize(myfunc, x0, args = (r,xiper,f1h0,f2h0),bounds=bounds)
print '0', res.x

res = minimize(myfunc, x0, args = (r,xipar,f1h1,f2h1),bounds=bounds)
print '1', res.x

x0 = res.x
res = minimize(myfunc, x0, args = (r,xi,f1h2,f2h2),bounds=bounds)
print '2', res.x

########################################################################################


#filename = '../perfil_%.2f_%.2f_%.2f_%.2f.%i' % (0.70,0.70,1.00,0.0,0)
#data = np.genfromtxt(filename)
#x = data[:,0]
#y = data[:,1]
#ax.plot(x,y/y[0]*1.4e4,'-',c='c')


#x,y = np.genfromtxt('../xilin_proyectada.dat',usecols=(0,1),unpack=True)
#ax.plot(x,y,'-g',label='xilin',)


#def integrando(x,r,alpha):
#    y = 1.0/np.power(np.sqrt(r*r + x*x),alpha)/np.power(1.0 + np.sqrt(r*r + x*x),3.0-alpha)
#    return y
#
#from scipy.integrate import quad
#n = 20
#rvec = np.linspace(-1,1,n)
#rvec = np.power(10.0,rvec)
#I = np.zeros(n)
#for r,i in zip(rvec,range(n)):
#    alpha = 2.0
#    I[i] = quad(integrando,0,np.inf,args=(r,alpha))[0]
#ax.plot(rvec,I/I[0]*1.4e4,'-',c='y')


ax1.set_xscale('log')
ax1.set_ylim(0.0,2.0)
ax1.set_xlabel(r'$r [h^{-1} Mpc]$',fontsize=18)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,23.0)
ax.set_ylim(6.0,5.0E3)
ax.set_ylabel(r'$\omega (r_{p})$',fontsize=18)

ax.set_title('Mass range: %s Magnitude: %2d'%(mass_range,mag))
rcParams['legend.frameon']=True
ax.legend()
ax.grid(True)

thismanager = plt.get_current_fig_manager()
#thismanager.set_window_title(suffix)
thismanager.window.setGeometry(0,0,800,853)
plt.show()
