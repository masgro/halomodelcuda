#!/usr/bin/env python

import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf
import os.path
import scipy.integrate as integrate
from scipy.interpolate import interp1d

def xilin(x):
      y = 0.804728
      y -= 0.963419*x
      y -= 0.408048*x*x
      y += 0.056507*x*x*x
      y += 0.096629*x*x*x*x
      y -= 0.049933*x*x*x*x*x
      y -= 0.032321*x*x*x*x*x*x
      return y


def xi(x,y):
          "Computa la funcion de correlacion"
          result = x/y - 1.0
          return result

def integrame(filename):
    from scipy.interpolate import interp1d

    file = filename
    if os.path.isfile(file):
        x,dd,dr  = np.genfromtxt(file,unpack='true',usecols=(0,1,2))
        y = dd/dr - 1.
    else:
        print 'NO EXISTE FILE: ',file
        return 1.0

    fileout = 'integrada_semianalitico.dat'
    fileout = open(fileout,"w")

    rmax = max(x)

    for rp in x:
        mask = (x >= rp)
        r = x[mask]

        if (len(r) < 2) :
            continue

        f = interp1d(r,y[mask])
        def myfunc(x) :
          y = f(np.sqrt(x*x+rp*rp))
          return(y)

        xxx = np.linspace(0.0,np.sqrt(rmax*rmax-rp*rp)-1.E-5,10000)
        ff = integrate.simps(myfunc(xxx),xxx)

        fileout.write("%f %f \n" % (rp,2.0*ff))
    fileout.close()



f = plt.figure()
axx = f.add_subplot(111)



def integra(folder,mmin,mmax,ang,bc,ab,alb,alc,mag,typ,t):

      if t == 0 :
        nlines = 100
      if t == 1 :
        nlines = 100

      ym1h = np.zeros((3, nlines))
      ym2h = np.zeros((3, nlines))
      ym   = np.zeros((3, nlines))
      xx   = np.zeros(nlines)

      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,bc,ab,alb,alc,mag)

      filename = folder+'/funcorr_1h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x1 = x[:,0]
        ym1h[0,:] = x[:,1]
        ym1h[1,:] = x[:,3]
        ym1h[2,:] = x[:,5]
        xx = x1
      else :
          print 'File not found: ',filename
      
      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,1.00,1.00,alb,alc,mag)
      filename = folder+'/funcorr_2h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x2 = x[:,0]
        ym2h[0,:] = x[:,1]
        ym2h[1,:] = x[:,3]
        ym2h[2,:] = x[:,5]
        xx = x2
      ###ym2h[0,:] = x[:,1]*(erf((np.log10(x2) - np.log10(r0))*3.0)*0.5 + 0.5)*factor*factor
      else :
          print 'File not found: ',filename
      
      ym = ym1h + ym2h


      xtmp = xx
      ytmp = ym[0,:]
      y1tmp = ym[1,:]
      y2tmp = ym[2,:]


      rmax = max(xx)
      #for rp in xtmp:
      #    mask = (xtmp > rp)
      #    r = xtmp[mask]
      #    xi = ytmp[mask]
      #    y1 = 2.0*integrate.simps(xi/np.sqrt(r*r-rp*rp),x=r)
      ##    y1 = integrate.trapz(xi/np.sqrt(r*r-rp*rp), x=r)
      #    print rp,y1


      filename = 'proyectada_integrada_'+suffix+'.dat'
      file = open(filename,"w")


      xxtmp = xtmp
      for rp in xxtmp:
          mask = (xtmp >= rp)
          r = xtmp[mask]

          if (len(r) < 2) :
              continue

          f = interp1d(r,ytmp[mask])
          def myfunc(x) :
            y = f(np.sqrt(x*x+rp*rp))
            return(y)

          xxx = np.linspace(0.0,np.sqrt(rmax*rmax-rp*rp)-1.E-5,10000)
          res0 = integrate.simps(myfunc(xxx),xxx)

          f = interp1d(r,y1tmp[mask])
          def myfunc1(x) :
            y = f(np.sqrt(x*x+rp*rp))
            return(y)
          res1 = integrate.simps(myfunc1(xxx),xxx)

          f = interp1d(r,y2tmp[mask])
          def myfunc2(x) :
            y = f(np.sqrt(x*x+rp*rp))
            return(y)
          res2 = integrate.simps(myfunc2(xxx),xxx)

          ##if t == 0 :
          ##    rmin = 10.0
          ##    rmax = 10.03
          ##if t == 1 :
          ##    rmin = 10.0
          ##    rmax = 10.3

          ##if(rp >= rmin) & (rp < rmax):
          ##  axx.plot(r,xi,'o')

          #factor = r/np.sqrt(r*r-rp*rp)
          #xi = ytmp[mask]*factor
          ##r1 = r[0:-1]
          ##r2 = r[1:]
          ##dr = r2-r1
          ##f = (xi[0:-1]+xi[1:])/2.0
          ##f = dr*f
          ##ff = sum(f)

          ##xi = y1tmp[mask]*factor
          ##f = (xi[0:-1]+xi[1:])/2.0
          ##f = dr*f
          ##ff1 = sum(f)

          ##xi = y2tmp[mask]*factor
          ##f = (xi[0:-1]+xi[1:])/2.0
          ##f = dr*f
          ##ff2 = sum(f)

          file.write("%f %f %f %f\n" % (rp,2.0*res0,2.0*res1,2.0*res2))

      file.close()

def plotfun(folder,mmin,mmax,ang,bc,ab,alb,alc,mag,typ,color,t):

      if t == 0 :
        nlines = 100
      if t == 1 :
        nlines = 100
      ym1h = np.zeros((3, nlines))
      ym2h = np.zeros((3, nlines))
      ym   = np.zeros((3, nlines))
      xx   = np.zeros(nlines)

      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,bc,ab,alb,alc,mag)

      filename = folder+'/funcorr_1h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x1 = x[:,0]
        ym1h[0,:] = x[:,1]
        ym1h[1,:] = x[:,3]
        ym1h[2,:] = x[:,5]
        xx = x1
      else :
          print 'NO existe file: ',filename
          ym1h = 0.0
       
      
      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,1.00,1.00,alb,alc,mag)
      filename = folder+'/funcorr_2h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x2 = x[:,0]
        ym2h[0,:] = x[:,1]
        ym2h[1,:] = x[:,3]
        ym2h[2,:] = x[:,5]
        xx = x2
        #ym2h[:,:] *= (erf((np.log10(x2) - np.log10(0.8))*3.0)*0.5 + 0.5)
      else :
          print 'NO existe file: ',filename
          ym2h = 0.0
      
      ym = ym1h + ym2h

      label = [r'$\xi_{a}$',r'$\xi_{b}$',r'$\xi_{c}$']
      colors = ['black','blue','green']
      for i in range(0,3):
        y = ym[i,:]
        ax1.plot(xx,y,'-',color=colors[i],label=label[i])
        #y = ym1h[i,:]
        #ax1.plot(xx,y, ':', color=colors[i])
        #y = ym2h[i,:]
        #ax1.plot(xx,y, ':', color=colors[i])
      
      filename = 'xi_'+suffix+'.dat'
      f = open(filename,'w')
      for i in range(len(xx)):
          f.write("%f %f\n" % (xx[i],ym[2,i]))
      f.close()



rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

f = plt.figure(figsize=(10,10))

ax1 = f.add_subplot(111)

ax1.set_xscale('log')
ax1.set_yscale('log')

c = ['b','r','g']
lt = ['k-.','k','k--']

#_bc = np.linspace(0.1,1.0,20)
#_ab = np.linspace(0.1,1.0,20)
#
#for bc in _bc[0:10]:
#  for ab in _ab[0:10]:
#    plotfun('../outputs/chi',13.50,14.00,45,bc,ab,1.00,1.00,19,'CG','red',0)

plotfun('../outputs',12.97,13.33,45,1.00,1.00,1.40,1.40,17,'CG','red',0)
#integra('../outputs/chi',13.74,15.00,45,0.80,0.80,1.40,1.40,17,'CG',0)

#_bc = np.linspace(0.1,1.4,50)
#_ab = np.linspace(0.1,1.4,50)
#
#for bc in _bc[0:10]:
#  for ab in _ab[0:10]:
#    plotfun(13.50,14.00,45,0.99,0.99,bc,ab,19,'CG','red')


###filename = '../1d/funcorr_1h_13.50-14.00_19.CG.dat'
###x1,ym1h = np.genfromtxt(filename,unpack=True,usecols=(0,1))
###filename = '../1d/funcorr_2h_13.50-14.00_19.CG.dat'
###x2,ym2h = np.genfromtxt(filename,unpack=True,usecols=(0,1))
###xx = x2
###ym = ym1h + ym2h
###ax1.plot(xx,ym,'-',lw=2,color='k',label='iso')
#####ax1.plot(xx,ym1h, ':',lw=1,color='k')
#####ax1.plot(xx,ym2h, ':',lw=1,color='k')


filename = 'Semianalitico/xibox.-17.12.97-13.33_mas11.dat'
x,dd,dr,dd1,dd2,dd3  = np.genfromtxt(filename,unpack='true',usecols=(0,1,2,3,4,5))

y = dd/dr - 1.
ax1.plot(x,y,'o',color='black',ms=5,lw=1,label=filename)

y = dd1/(dr*0.3) - 1.
ax1.plot(x,y,'o',color='red',ms=5,lw=1)
y = dd2/(dr*0.3) - 1.
ax1.plot(x,y,'o',color='blue',ms=5,lw=1)
y = dd3/(dr*0.3) - 1.
ax1.plot(x,y,'o',color='green',ms=5,lw=1)
#
integrame(filename)

ax1.set_xlim(0.1,50.0)
ax1.set_ylim(1.0E-2,1.0E4)
ax1.set_ylabel(r'$\xi (r)$')
ax1.set_xlabel(r'$r [h^{-1} Mpc]$')

rcParams['legend.frameon']=False
ax1.legend()
#ax1.grid(True)

plt.show()
