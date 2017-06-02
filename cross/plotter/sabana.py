from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.special import erf
from scipy import interpolate
from matplotlib import cm
from numpy.random import uniform, seed
from matplotlib.mlab import griddata
import math
import matplotlib.animation as animation
from scipy import integrate

C = 0.0
D = 0.0

def g(x,y):
    global C,D
    f  = np.exp(-x*x*0.5/C/C)
    f  = f + D*np.exp(-np.power(x-math.pi*0.5,2)*0.5/C/C)*np.exp(-y*y*0.5/C/C)
    return f

rc('text', usetex=False)
rcParams.update({'font.size': 12})

fig = plt.figure(figsize=(10,8))
ax = fig.gca(projection='3d')

ims = []
for i in range(20):
    C = 0.7
    D = 0.1 + float(i)*0.1
    
    area = integrate.nquad(g, [[0.0,math.pi*0.5],[0.0,math.pi*0.5]])[0]

    X = np.linspace(0.0,math.pi*0.5,200)
    Y = np.linspace(0.0,math.pi*0.5,200)
    X, Y = np.meshgrid(X, Y)
    Z = g(X,Y)/area
    
    surf = ax.plot_surface(X, Y, Z, rstride=4, cstride=4, alpha=1.0, cmap='hsv', linewidth=0.1)
    
    #cset = ax.contour(X, Y, Z, zdir='z', offset=0)

    ims.append([surf])


ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True,repeat_delay=100)


ax.zaxis.set_major_locator(LinearLocator(5))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel(r'$\theta$')
ax.set_xlim(0, math.pi*0.5)
ax.set_ylabel(r'$\phi$')
ax.set_ylim(0, math.pi*0.5)
ax.set_zlabel(r'$f(\vec{y}_{1,2},\vec{\epsilon}_1)$')
ax.set_zlim(0, 1)

ani.save("dynamic_images.mp4")

plt.show()
