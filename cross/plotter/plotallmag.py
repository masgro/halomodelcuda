import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf
import os.path

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

def plotfun(folder,mmin,mmax,ang,bc,ab,alb,alc,mag,typ,color,ax,lw,factor):
      nlines = 50
      ym1h = np.zeros((3, nlines))
      ym2h = np.zeros((3, nlines))
      ym   = np.zeros((3, nlines))
      xx   = np.zeros(nlines)

      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,bc,ab,1.0,1.0,mag)

      filename = folder+'/funcorr_1h_'+suffix+'.'+typ+'.dat'
      print filename
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x1 = x[:,0]
        ym1h[0,:] = x[:,1]#*factor
        ym1h[1,:] = x[:,3]#*factor
        ym1h[2,:] = x[:,5]#*factor
        xx = x1

      suffix = '%.2f-%.2f_%02d_%.2f_%.2f_%.2f_%.2f_%2d' % (mmin,mmax,ang,0.99,0.99,alb,alc,mag)
      
      filename = folder+'/funcorr_2h_'+suffix+'.'+typ+'.dat'
      print filename
      if os.path.isfile(filename):
        x = np.genfromtxt(filename)
        x2 = x[:,0]
        ym2h[0,:] = x[:,1]*factor
        ym2h[1,:] = x[:,3]*factor
        ym2h[2,:] = x[:,5]*factor
        xx = x2
        #ym2h *= (erf((np.log10(x2) - np.log10(1.0))*3.0)*0.5 + 0.5)
      
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

      r0 = 6.
      g =  1.5
      rr=np.linspace(0.01,100,1000)
      xi = np.power(r0/rr,g)
      ax.plot(rr,xi,'-ok',lw=2)


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

  
def plot_xi(filename,ax):
    ###x = np.genfromtxt(filename)
    ###x1 = x[:,0]/1000.0
    ###y4 = xi(x[:,5],x[:,6])
    ###ax.plot(x1,y4,marker='o')
          
    if os.path.isfile(filename):
        x,dd,dr,dd1,dd2,dd3  = np.genfromtxt(filename,unpack='true',usecols=(0,1,2,3,4,5))
        y = dd/dr - 1.
        ax.plot(x,y,'o',color='black',ms=5,lw=1,label='19')
        y = dd1/(dr*0.3) - 1.
        ax.plot(x,y,'o',color='red',ms=5,lw=1,label='19')
        y = dd2/(dr*0.3) - 1.
        ax.plot(x,y,'o',color='blue',ms=5,lw=1,label='19')
        y = dd3/(dr*0.3) - 1.
        ax.plot(x,y,'o',color='green',ms=5,lw=1,label='19')


rc('text', usetex=False)
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


f = plt.figure(figsize=(10,10))

gs = gridspec.GridSpec(3,2)
gs.update(left=0.1, right=0.96, hspace=0,wspace=0)

ax = []
ax.append(plt.subplot(gs[0,0]))
ax.append(plt.subplot(gs[0,1]))
ax.append(plt.subplot(gs[1,0]))
ax.append(plt.subplot(gs[1,1]))
ax.append(plt.subplot(gs[2,0]))
ax.append(plt.subplot(gs[2,1]))


m1 = 13.0
m2 = 13.5
magv = [16,17,18,19,20,21]
bcv = [0.76,0.69,0.60,0.56,0.49,0.39]
abv = [0.71,0.69,0.67,0.65,0.63,0.63]
ala = [0.76,0.82,0.71,0.82,0.76,0.71]
alb = [0.26,0.34,0.42,0.29,0.29,0.34]

ni = 6
for i in range(ni):
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    ax[i].set_xlim(0.1,40.0)
    ax[i].set_ylim(5.0E-2,1.0E4)
    label = 'Mag %.02f' % magv[i]
    ax[i].text(7,3000, label)


ax[0].set_ylabel(r'$\xi (r)$')
ax[2].set_ylabel(r'$\xi (r)$')
ax[2].set_xlabel(r'$r [h^{-1} Mpc]$')
ax[3].set_xlabel(r'$r [h^{-1} Mpc]$')

xticklabels = ax[0].get_xticklabels()+ax[1].get_xticklabels()
setp(xticklabels, visible=False)
yticklabels = ax[1].get_yticklabels()+ax[3].get_yticklabels()
setp(yticklabels, visible=False)

f.suptitle('%.02f - %.02f' % (m1,m2))

for i in range(ni):
    plotfun('../outputs/chi',m1,m2,45,bcv[i],abv[i],ala[i],alb[i],magv[i],'CG','red',ax[i],3,1.0)
    #plotfun('../outputs/chi',m1,m2,90,0.99,0.99,0.90,0.30,16,'HM','red',ax[i],1,1.0)

    filename = 'Semianalitico/xibox.-'+ '%2d' % magv[i] +'.'+'%.02f' % m1+'-'+'%.02f' % m2+'.mas11.dat'
    plot_xi(filename,ax[i])


plt.show()
