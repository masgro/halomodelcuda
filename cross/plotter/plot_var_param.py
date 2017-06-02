import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf
from scipy import interpolate

def plot_hm_tipo(bc,ab,al_b,al_c,tipo):
          sbc = "%.2f" % bc
          sab = "%.2f" % ab
          sal_b = "%.2f" % al_b
          sal_c = "%.2f" % al_c
          file = 'funcorr_1h_12.00-16.00_01_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'.dat'
          data = np.genfromtxt(file)
          x = data[:,0]
          y1h = data[:,tipo]
          if tipo == 1:
            t = 'A'
          if tipo == 2:
            t = 'B'
          if tipo == 3:
            t = 'C'
          file = 'funcorr_2h_12.00-16.00_01_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'.dat.'+t
          data = np.genfromtxt(file)
          x = data[:,0]
          y2h = data[:,tipo]
          return x,y1h+y2h

rc('text', usetex=False)
rcParams.update({'font.size': 8})

f = plt.figure(figsize=(8,4))

gs = gridspec.GridSpec(2,4)
gs.update(left=0.1, right=0.96, hspace=0.00, wspace=0.00)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])
ax3 = plt.subplot(gs[3])
ax4 = plt.subplot(gs[4])
ax5 = plt.subplot(gs[5])
ax6 = plt.subplot(gs[6])
ax7 = plt.subplot(gs[7])

ax0.set_xscale('log')
ax0.set_yscale('log')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax4.set_xscale('log')
ax4.set_yscale('log')
ax5.set_xscale('log')
ax5.set_yscale('log')
ax6.set_xscale('log')
ax6.set_yscale('log')
ax7.set_xscale('log')
ax7.set_yscale('log')

i = 0

m1 = ['12.00','12.50','13.00','13.50','14.00','14.50']
m2 = ['12.50','13.00','13.50','14.00','14.50','15.00']

xbest = [0.85, 0.00, 0.00, 0.00, 0.00, 0.57] 
ybest = [0.45, 0.00, 0.00, 0.00, 0.00, 0.33]

m1 = m1[i]
m2 = m2[i]
xbest = xbest[i]
ybest = ybest[i]

angulo = '01'
term   = '1h'

ax = ax0
A = 0.80
B = 0.80
C = 0.91
D = 0.35
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{I}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)



ax = ax1
A = 0.80
B = 0.20
C = 0.91
D = 0.35
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{II}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)



ax = ax2
A = 0.20
B = 0.80
C = 0.91
D = 0.35
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{III}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)



ax = ax3
A = 0.99
B = 0.99
C = 0.91
D = 0.35
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{IV}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)




ax = ax4
A = 0.80
B = 0.80
C = 0.91
D = 0.10
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{V}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)



ax = ax5
A = 0.80
B = 0.80
C = 0.91
D = 1.00
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{VI}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)



ax = ax6
A = 0.80
B = 0.80
C = 0.50
D = 0.35
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{VII}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)




ax = ax7
A = 0.80
B = 0.80
C = 3.00
D = 0.35
ax.text(3.0, np.power(10.0,3.0), 'A = ' + "%.2f " % A)
ax.text(3.0, np.power(10.0,2.6), 'B = ' + "%.2f " % B)
ax.text(3.0, np.power(10.0,2.2), 'C = ' + "%.2f " % C)
ax.text(3.0, np.power(10.0,1.8), 'D = ' + "%.2f " % D)
ax.text(0.15, 0.20, r'$\mathrm{VIII}$')
x,y = plot_hm_tipo(A,B,C,D,1)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '-', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,2)
tck = interpolate.splrep(x,y,s=1)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, ':', color='black',linewidth = 1)

x,y = plot_hm_tipo(A,B,C,D,3)
tck = interpolate.splrep(x,y,s=0)
y = interpolate.splev(x,tck,der=0)
ax.plot(x, y, '--', color='black',linewidth = 1)




majorFormatter = FormatStrFormatter('%.1f')

ax0.set_xlim(0.1,30.00)
ax0.set_ylim(5.0E-2,5.0E3)
ax0.set_ylabel(r'$\xi (r)$')

ax1.set_xlim(0.1,30.00)
ax1.set_ylim(5.0E-2,5.0E3)

ax2.set_xlim(0.1,30.00)
ax2.set_ylim(5.0E-2,5.0E3)

ax3.set_xlim(0.1,30.00)
ax3.set_ylim(5.0E-2,5.0E3)

ax4.set_xlim(0.1,30.00)
ax4.set_ylim(5.0E-2,5.0E3)
ax4.set_xlabel(r'$r [h^{-1} Mpc]$')
ax4.set_ylabel(r'$\xi (r)$')
ax4.xaxis.set_major_formatter(majorFormatter) 

ax5.set_xlim(0.1,30.00)
ax5.set_ylim(5.0E-2,5.0E3)
ax5.set_xlabel(r'$r [h^{-1} Mpc]$')
ax5.xaxis.set_major_formatter(majorFormatter) 

ax6.set_xlim(0.1,30.00)
ax6.set_ylim(5.0E-2,5.0E3)
ax6.set_xlabel(r'$r [h^{-1} Mpc]$')
ax6.xaxis.set_major_formatter(majorFormatter) 

ax7.set_xlim(0.1,30.00)
ax7.set_ylim(5.0E-2,5.0E3)
ax7.set_xlabel(r'$r [h^{-1} Mpc]$')
ax7.xaxis.set_major_formatter(majorFormatter) 

xticklabels = ax0.get_xticklabels()+ax1.get_xticklabels()+ax2.get_xticklabels()+ax3.get_xticklabels()
setp(xticklabels, visible=False)
yticklabels = ax1.get_yticklabels()+ax2.get_yticklabels()+ax3.get_yticklabels()
setp(yticklabels, visible=False)
yticklabels = ax5.get_yticklabels()+ax6.get_yticklabels()+ax7.get_yticklabels()
setp(yticklabels, visible=False)

plt.show()

filename = 'variacion_de_parametros.ps'
savefig(filename)
