import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf
import os.path

from scipy import integrate

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


def integrame(folder,mmin,mmax,ang,bc,ab,alb,alc,mag,typ):
      nlines = 50
      ym1h = np.zeros((3, nlines))
      ym2h = np.zeros((3, nlines))
      ym   = np.zeros((3, nlines))
      xx   = np.zeros(nlines)

      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,bc,ab,alb,alc,mag)

      filename = folder+'/funcorr_1h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        print filename
        x = np.genfromtxt(filename)
        x1 = x[:,0]
        ym1h[0,:] = x[:,1]
        ym1h[1,:] = x[:,3]
        ym1h[2,:] = x[:,5]
        xx = x1
      
      filename = folder+'/funcorr_2h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        print filename
        x = np.genfromtxt(filename)
        x2 = x[:,0]
        ym2h[0,:] = x[:,1]
        ym2h[1,:] = x[:,3]
        ym2h[2,:] = x[:,5]
        xx = x2
      ###ym2h[0,:] = x[:,1]*(erf((np.log10(x2) - np.log10(r0))*3.0)*0.5 + 0.5)*factor*factor
      
      ym = ym1h + ym2h


      xtmp = xx
      ytmp = ym[0,:]

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
          mask = (xtmp > rp)
          r = xtmp[mask]
          xi = ytmp[mask]
          xi /= np.sqrt(r*r-rp*rp)
          xi *= r


          r1 = r[0:-1]
          r2 = r[1:]
          dr = r2-r1
          f = (xi[0:-1]+xi[1:])/2.0
          f = dr*f
          ff = sum(f)

          file.write("%f %f\n" % (rp,2.0*ff))

      file.close()



def integrame_semianalitico(folder,filename):
    from scipy.interpolate import interp1d


    file = folder+filename
    if os.path.isfile(file):
        x,dd,dr  = np.genfromtxt(file,unpack='true',usecols=(0,1,2))
        y = dd/dr - 1.
    else:
        print 'NO EXISTE FILE: ',file
        return 1.0

    print "Semianalitico"
    fileout = 'integrado_'+ filename
    fileout = open(fileout,"w")

    xtmp = x
    for rp in xtmp:
        mask = (x > rp)
        r = x[mask]
        xi = y[mask]
        xi /= np.sqrt(r*r-rp*rp)
        xi *= r

        r1 = r[0:-1]
        r2 = r[1:]
        dr = r2-r1
        f = (xi[0:-1]+xi[1:])/2.0
        f = dr*f
        ff = sum(f)
        fileout.write("%f %f \n" % (rp,2.0*ff))
    fileout.close()


def plotfun(folder,mmin,mmax,ang,bc,ab,alb,alc,mag,typ,color,ax,lw,factor):
      nlines = 50
      ym1h = np.zeros((3, nlines))
      ym2h = np.zeros((3, nlines))
      ym   = np.zeros((3, nlines))
      xx   = np.zeros(nlines)

      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,bc,ab,alb,alc,mag)

      filename = folder+'/funcorr_1h_'+suffix+'.'+typ+'.dat'
      print filename
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x1 = x[:,0]
        ym1h[0,:] = x[:,1]#*factor
        ym1h[1,:] = x[:,3]#*factor
        ym1h[2,:] = x[:,5]#*factor
        xx = x1
      
      filename = folder+'/funcorr_2h_'+suffix+'.'+typ+'.dat'
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x2 = x[:,0]
        ym2h[0,:] = x[:,1]*factor
        ym2h[1,:] = x[:,3]*factor
        ym2h[2,:] = x[:,5]*factor
        xx = x2
      ###ym2h[0,:] = x[:,1]*(erf((np.log10(x2) - np.log10(r0))*3.0)*0.5 + 0.5)*factor*factor
      
      ym = ym1h + ym2h

      label = [r'$\xi_{a}$',r'$\xi_{b}$',r'$\xi_{c}$']
      colors = ['red','blue','green']
      for i in range(0,3):
        y = ym[i,:]
        ax.plot(xx,y,'-',color=colors[i],label=label[i],lw=lw)
        y = ym1h[i,:]
        ax.plot(xx,y, ':', color=colors[i])
        y = ym2h[i,:]
        ax.plot(xx,y, ':', color=colors[i])

      ax.plot(xx,np.power(xx/8.0,-2.1),'--')


def plot_hm(mmin,mmax,r0):
          c = ['b','r','g']
          
          nfiles = 5
          nlines = 100
          name = 'cono'
          angulo = '90'

          ym = np.zeros((3, nfiles, nlines))
          y1h = np.zeros((3, nfiles, nlines))
          y2h = np.zeros((3, nfiles, nlines))
          
          factor = 0.9
          norma_merchan = 3.579438e-02
          
          b = 0
          while b < nfiles:
            a = 0
            while a < 3:
              file = 'data/funcorr_'+str(a)+'_1h_'+mmin+'-'+mmax+'_'+name+'_'+angulo+'_'+str(b).zfill(2)+'.dat'
              x = np.genfromtxt(file)
              x1 = x[:,0]
              y1 = x[:,1]*factor/norma_merchan
            
              file = 'data/funcorr_'+str(a)+'_2h_'+mmin+'-'+mmax+'_'+name+'_'+angulo+'_'+str(b).zfill(2)+'.dat'
              x = np.genfromtxt(file)
              x2 = x[:,0]
              y2 = x[:,1]*factor*factor*(erf((np.log10(x2) - np.log10(r0))*3.0)*0.5 + 0.5)*0.7/norma_merchan/norma_merchan
          
              ym[a,b,:] = (y1 + y2)
              y1h[a,b,:] = y1
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
            ax.plot(x1,y, lt[a], color='black',linewidth = 1)

            y = y1h[a,:,:]
            y.reshape(nfiles,nlines)
            y = np.median(y, axis=0)
            yy[a,:] = y
            ax.plot(x1,y, ':', color='black',linewidth = 1)
          
            y = y2h[a,:,:]
            y.reshape(nfiles,nlines)
            y = np.median(y, axis=0)
            yy[a,:] = y
            ax.plot(x1,y, ':', color='black',linewidth = 1)
          
            a += 1

  
def plot_xi(folder,filename,ax):
    ###x = np.genfromtxt(filename)
    ###x1 = x[:,0]/1000.0
    ###y4 = xi(x[:,5],x[:,6])
    ###ax.plot(x1,y4,marker='o')
    
    filename = folder+filename

    if os.path.isfile(filename):
        x,dd,dr,dd1,dd2,dd3  = np.genfromtxt(filename,unpack='true',usecols=(0,1,2,3,4,5))
        y = dd/dr - 1.
        ax.plot(x,y,'o',color='black',ms=5,lw=1,label='19')
        #y = dd1/(dr*0.3) - 1.
        #ax.plot(x,y,'o',color='red',ms=5,lw=1,label='19')
        #y = dd2/(dr*0.3) - 1.
        #ax.plot(x,y,'o',color='blue',ms=5,lw=1,label='19')
        #y = dd3/(dr*0.3) - 1.
        #ax.plot(x,y,'o',color='green',ms=5,lw=1,label='19')


rc('text', usetex=False)
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


f = plt.figure(figsize=(10,10))

gs = gridspec.GridSpec(2,2)
gs.update(left=0.1, right=0.96, hspace=0,wspace=0)

ax = []
ax.append(plt.subplot(gs[0,0]))
ax.append(plt.subplot(gs[0,1]))
ax.append(plt.subplot(gs[1,0]))
ax.append(plt.subplot(gs[1,1]))

ni = 4
mmin = 11.50
mmax = 15.50
dm = (mmax - mmin)/ni
#m1v = [12.50,13.00,13.50,14.00]
#m2v = [13.00,13.50,14.00,14.50]


m1v = [12.50, 12.78, 13.12, 13.55]
m2v = [12.78, 13.12, 13.55, 15.00]

for i in range(ni):
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    ax[i].set_xlim(0.1,40.0)
    ax[i].set_ylim(5.0E-2,1.0E4)
    label = '%.02f-%.02f' % (m1v[i],m2v[i])
    ax[i].text(7,3000, label)


ax[0].set_ylabel(r'$\xi (r)$')
ax[2].set_ylabel(r'$\xi (r)$')
ax[2].set_xlabel(r'$r [h^{-1} Mpc]$')
ax[3].set_xlabel(r'$r [h^{-1} Mpc]$')

xticklabels = ax[0].get_xticklabels()+ax[1].get_xticklabels()
setp(xticklabels, visible=False)
yticklabels = ax[1].get_yticklabels()+ax[3].get_yticklabels()
setp(yticklabels, visible=False)

integrame('../outputs/chi',13.00,13.50,90,0.99,0.99,1.40,1.40,19,'CG')
#

#filename = 'Semianalitico/xibox.-19.'+'%.02f' % 13.0 +'-'+'%.02f' % 13.5 + '.mas11.dat'
#integrame_semianalitico(filename)

mag = 21

for i in range(ni):
    m1 = m1v[i]
    m2 = m2v[i]
    plotfun('../outputs/chi',m1,m2,90,0.99,0.99,1.40,1.40,mag,'CG','red',ax[i],3,1.0)
    plotfun('../outputs/chi',m1,m2,90,0.99,0.99,1.40,1.40,mag,'HM','red',ax[i],1,1.0)
    integrame('../outputs/chi',m1,m2,90,0.99,0.99,1.40,1.40,mag,'CG')

    folder = 'Semianalitico/'
    filename = 'xibox.-21.'+'%.02f' % m1+'-'+'%.02f' % m2+'.dat'
    plot_xi(folder,filename,ax[i])

    integrame_semianalitico(folder,filename)




plt.show()
