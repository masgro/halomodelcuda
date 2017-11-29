#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import math as m
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.integrate import quad

def xilin(x):
      y  = 0.778235
      y -= 0.961176*x
      y -= 0.344563*x*x
      y += 0.040576*x*x*x
      y += 0.082509*x*x*x*x
      y -= 0.049040*x*x*x*x*x
      y -= 0.031135*x*x*x*x*x*x
      return y


def integrand(x,r):
      return x*np.power(10.0,xilin(np.log10(x)))/np.sqrt(x*x + r*r)


xmin = -2.0
xmax = +2.0
ntab = 100
x = np.linspace(xmin,xmax,num=ntab)
y = np.arange(ntab, dtype=np.float)

for i in np.arange(ntab): 
  r = np.power(10.0,x[i])
  y[i] = quad(integrand, r, 10000.0, args=(r))[0]


rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(111)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.01,100.0)
ax1.set_ylim(1.0E-2,1.0E2)

ax1.plot(np.power(10.0,x),np.power(10.0,xilin(x)),'black')
ax1.plot(np.power(10.0,x),y,'red')

z = np.polyfit(x, np.log10(y), 6)
print z

yy  = z[0]*x*x*x*x*x*x
yy += z[1]*x*x*x*x*x
yy += z[2]*x*x*x*x
yy += z[3]*x*x*x
yy += z[4]*x*x
yy += z[5]*x
yy += z[6] 

ax1.plot(np.power(10.0,x),np.power(10.0,yy),'blue')



plt.show()



