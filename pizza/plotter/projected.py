#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import math as m
from matplotlib import rcParams

def rho(r):
      y = 1./r/(1. + r)/(1.0 + r)
      return y

def f(r):
      y = np.zeros_like(r)

      indx = np.where(r < 1.0);
      r1 = r[indx]
      y[indx] = np.sqrt(1.0 - r1)/np.sqrt(1.0 + r1)
      y[indx] = -1.0 + 2.0/np.sqrt(1.0 - r1*r1)*np.arctanh(y[indx])
      y[indx] = y[indx]/(1.0 - r1*r1)

      indx = np.where(r > 1.0);
      r2 = r[indx]
      y[indx] = np.sqrt(r2 - 1.0)/np.sqrt(r2 + 1.0)
      y[indx] = 1.0 - 2.0/np.sqrt(r2*r2 - 1.0)*np.arctan(y[indx])
      y[indx] = y[indx]/(r2*r2 - 1.0)

      return y

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True
fig = plt.figure(figsize=(8,8))
ax1 = fig.add_subplot(111)
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_xlim(0.1,100.0)
#ax1.set_ylim(1.0E-6,1.0E4)


theta = m.pi/2.0
phi   = m.pi/2.0

phii = np.linspace(0.001,m.pi/2.0,10)
thetai = np.linspace(0.001,m.pi/2.0,10)
for i in range(len(phii)):
  for j in range(len(thetai)):
     theta = thetai[j]
     phi = phii[i]

     a = 0.25
     b = 0.5
     c = 1.0
     
     st = np.sin(theta)
     ct = np.cos(theta)
     sp = np.sin(phi)
     cp = np.cos(phi)
     
     st2 = st*st
     ct2 = ct*ct
     sp2 = sp*sp
     cp2 = cp*cp
     
     cc = c*c
     bb = b*b
     aa = a*a
     
     ff = st2*(cc/aa*cp2 + cc/bb*sp2) + ct2
     
     AM = ct2*(cc/aa*sp2 + cc/bb*cp2) + cc/aa*cc/bb*st2
     BM = ct*np.sin(2.0*phi)*(cc/aa - cc/bb)
     CM = cc/bb*sp2 + cc/aa*cp2
     
     psi = 0.5*np.arctan(BM/(AM-CM))
     
     SQRTABC = np.sqrt((AM-CM)*(AM-CM)+BM*BM)
     qx2 = 2.0*ff/(AM + CM - SQRTABC)
     qy2 = 2.0*ff/(AM + CM + SQRTABC)
     
     dc = 2.5044467e+04
     rhocrit = 2.7755001e+11
     
     x = np.linspace(-1.0, 2.0)
     x = np.power(10.0,x)
     y = np.zeros_like(x)
     z = y
     r = c*c*((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c))
     r = np.sqrt(r)
     #ax1.plot(x,dc*rhocrit/1.E13*rho(r),'black')
     
     y = x
     x = x*0.0
     z = x
     r = c*c*((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c))
     r = np.sqrt(r)
     #ax1.plot(y,dc*rhocrit/1.E13*rho(r),'red')
     
     z = y
     y = y*0.0
     x = y
     r = c*c*((x*x)/(a*a) + (y*y)/(b*b) + (z*z)/(c*c))
     r = np.sqrt(r)
     #ax1.plot(z,dc*rhocrit/1.E13*rho(r),'blue')
     
     x = np.linspace(-1.0, 2.0)
     x = np.power(10.0,x)
     
     r = x*x/qx2
     r = np.sqrt(r)
     #ax1.plot(x,2.0*dc*rhocrit*f(r)/np.sqrt(ff)/1.e13,'green')
     
     r = x*x/qy2
     r = np.sqrt(r)
     #ax1.plot(x,2.0*dc*rhocrit*f(r)/np.sqrt(ff)/1.e13,'orange')

abmedio = 1.00
bcmedio = 0.70
costheta = 1.0
phi = 0.0

colors = ['c','k','b','r']
bcv = [1.00,1.0,0.7,0.7]
abv = [1.00,0.7,0.7,0.5]
costhetav = [1.0,1.0,1.0,1.0]
phiv = [0.0,0.0,0.0,0.0]
for bcmedio,abmedio,c,costheta,phi in zip(bcv,abv,colors,costhetav,phiv):
    for i,style in zip(range(3),['-^','-o','--']):
        filename = '../perfil_%.2f_%.2f_%.2f_%.2f.%i' % (bcmedio,abmedio,costheta,phi,i)
        data = np.genfromtxt(filename)
        x = data[:,0]
        y = data[:,1]
        ax1.plot(x,y*1E4,style,c=c)

plt.ylim(1.E-11,1.e-9)
plt.show()
