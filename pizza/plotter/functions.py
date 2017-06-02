import numpy as np
import os.path
from scipy.special import erf,gamma

def plot_halo_model(folder,mag,mass_range,bc,ab,al_b,al_c,ax,ax1):
    nlines = 50
    ndir = 3
    ym1h = np.zeros((ndir, nlines))
    er1h = np.zeros((ndir, nlines))
    ym2h = np.zeros((ndir, nlines))
    er2h = np.zeros((ndir, nlines))
    ym   = np.zeros((ndir, nlines))
    er   = np.zeros((ndir, nlines))
    
    angulo = '00'
    sbc = "%.2f" % bc
    sab = "%.2f" % ab
    sal_b = "%.2f" % 1.00
    sal_c = "%.2f" % 1.00
    
    suffix = mass_range+'_'+angulo+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+'%2d'%mag+'.CG.dat'
    filename1h = folder+'funcorr_1h_'+suffix
    print filename1h
    
    sbc = "%.2f" % 1.00
    sab = "%.2f" % 1.00
    sal_b = "%.2f" % al_b
    sal_c = "%.2f" % al_c
    
    suffix = mass_range+'_'+angulo+'_'+sbc+'_'+sab+'_'+sal_b+'_'+sal_c+'_'+'%2d'%mag+'.CG.dat'
    filename2h = folder+'funcorr_2h_'+suffix
    print filename2h
    
    if os.path.isfile(filename1h):
        data = np.genfromtxt(filename1h,unpack=True)
        x = data[0]
        ym1h[0,:] = data[1]
        er1h[0,:] = data[2]
        ym1h[1,:] = data[3]
        er1h[1,:] = data[4]
        ym1h[2,:] = data[5]
        er1h[2,:] = data[6]
    else :
        print 'Does not exists file: ',filename1h
        return 0

    
    if os.path.isfile(filename2h):
        data = np.genfromtxt(filename2h,unpack=True)
        x = data[0]
        ym2h[0,:] = data[1]
        er2h[0,:] = data[2]
        ym2h[1,:] = data[3]
        er2h[1,:] = data[4]
        ym2h[2,:] = data[5]
        er2h[2,:] = data[6]
    else :
        print 'Does not exists file: ',filename2h
        return 0
    
    r0 = 0.9
    ym2h = (erf((np.log10(x) - np.log10(r0))*3.0)*0.5 + 0.5)*ym2h
    
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
    ax.plot(x,y,'-',lw=2,color=color)
    #ax.fill_between(x,y-error,y+error,alpha=0.5,edgecolor=color,facecolor=color)
    
    y = ym1h[0,:]
    ax.plot(x,y, ':', color='black')
    y = ym1h[1,:]
    ax.plot(x,y, ':', color='black')
    y = ym1h[2,:]
    ax.plot(x,y, ':', color='black')
    
    y = ym2h[0,:]
    ax.plot(x,y, ':', color='black')
    y = ym2h[1,:]
    ax.plot(x,y, ':', color='black')
    y = ym2h[2,:]
    ax.plot(x,y, ':', color='black')
    
    color = 'red'
    ax1.plot(x,ym[0,:]/ym[2,:], '-', lw=2,color=color)
    color = 'blue'
    ax1.plot(x,ym[1,:]/ym[2,:], '-', lw=2,color=color)
    color = 'black'
    ax1.plot(x,ym[2,:]/ym[2,:],'-',lw=2,color=color)
