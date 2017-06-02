import numpy as np
import math

def cvir(lgm):
  q = 1.020 - 0.109*(lgm - 12.0)
  q = np.power(10.0,q)
  return q

def rvir(lgm):
  q = np.log10(4.0/3.0*math.pi*deltavir*rhocrit)
  q = (lgm - q)/3.0
  q = np.power(10.0,q)
  return q

def deltae(ab,bc):
  q = 5.0*deltavir
  q = q*np.power(1.0/(bc*bc*ab),0.75)
  return q

def ce(cvir):
  q = cvir*0.5
  return q

def dc(deltae, ce):
  q = deltae*ce*ce*ce/3.0
  q = q/(np.log(1.0 + ce) - ce/(1 + ce))
  return q

def factores(theta,phi,ab,bc):
  ct = np.cos(theta)
  st = np.sin(theta)
  cp = np.cos(phi)
  sp = np.sin(phi)

  st2 = st*st
  ct2 = ct*ct
  sp2 = sp*sp
  cp2 = cp*cp

  f = st2*(cp2/(ab*ab) + sp2)/(bc*bc) + ct2

  A = ct2*(sp2/(ab*ab) + cp2)/(bc*bc) + st2/(bc*bc)/(bc*bc)/(ab*ab)
  B = ct*np.sin(2.0*phi)*(1.0/(ab*ab) - 1.0)/(bc*bc)
  C = (sp2 + cp2/(ab*ab))/(bc*bc)

  qx2 = 2.0*f/(A + C - np.sqrt(np.power(A-C,2.0) + B*B))
  qy2 = 2.0*f/(A + C + np.sqrt(np.power(A-C,2.0) + B*B))

  return f,A,B,C,qx2,qy2


def f_projected(r):

  q = r
  rr = 1.0 - r*r
  i = np.where(rr > 0.0)
  if len(i) > 0:
    q[i] = np.sqrt((1.0 - r[i])/(1.0 + r[i]))
    q[i] = -1.0 + 2.0/np.sqrt(rr[i])*np.arctanh(q[i])
    q[i] = q[i]/(rr[i])

  i = np.where(rr < 0.0)
  if len(i) > 0:
    q[i] = np.sqrt((r[i] - 1.0)/(r[i] + 1.0))
    q[i] = 1.0 - 2.0/np.sqrt(-rr[i])*np.arctan(q[i])
    q[i] = q[i]/(-rr[i])

  i = np.where(rr == 0.0)
  if len(i) > 0:
    q[i] = 1.0/3.0

  return q



def projected_profile(x,lgm,ab,bc,theta,phi,term):
  cv = cvir(lgm)
  rv = rvir(lgm)
  d = deltae(ab,bc)
  c = ce(cv)

  rs = rv/cv

  f,A,B,C,qx2,qy2 = factores(theta,phi,ab,bc)

  xx = x[0,:]
  yy = x[1,:]

  if term==1:
    r = xx*xx/qx2 + yy*yy/qy2
  elif term==2:
    r = (A*xx*xx+B*xx*yy+C*yy*yy)/f

  r = np.sqrt(r)
  r = r/rs

  rho = 2.0*dc(d,c)*rhocrit
  rho = rho/np.sqrt(f)
  rho = rho*f_projected(r)
  return rho

def profile(x,lgm,ab,bc):
  cv = cvir(lgm)
  rv = rvir(lgm)
  d = deltae(ab,bc)
  c = ce(cv)

  rs = rv/cv

  xx = x[0,:]
  yy = x[1,:]
  zz = x[2,:]

  r = zz*zz + yy*yy/bc/bc + xx*xx/(ab*ab*bc*bc);
  r = np.sqrt(r)
  r = r/rs

  rho = dc(d,c)*rhocrit
  rho = rho/r
  rho = rho/np.power(1.0+r,2)
  return rho

import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf
import matplotlib.cm as cm


rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

f = plt.figure(figsize=(8,8))

ax1 = f.add_subplot(111)

ax1.set_xscale('log')
ax1.set_yscale('log')
deltavir=95.4
rhomedio=7.1608e+10
rhocrit=2.7755e11

amax = math.pi*0.5
angles = np.arange(0.0,amax,amax/20.0)
r = np.arange(-1.0,2.0,0.1)
r = np.power(10.0,r)

for angle in angles:
  psi = 0.0
  x1 = r*np.cos(psi)
  x2 = r*np.sin(psi)

  x = np.array([x1,x2],np.float32)

  theta = math.pi/4.0
  phi   = angle
  lgm   = 13.0
  ab    = 1.0
  bc    = 0.5

  y1 = projected_profile(x,lgm,ab,bc,theta,phi,1)/np.power(10.0,lgm)
  color = angle/amax
  ax1.plot(r,y1, '-', c=(color,0,0,1.0),lw=2)

  y2 = projected_profile(x,lgm,ab,bc,theta,phi,2)/np.power(10.0,lgm)
  color = angle/amax
  ax1.plot(r,y2, '-', c=(0,0,color,1.0),lw=2)


  x1 = r
  x2 = r*0.0
  x3 = r*0.0
  x = np.array([x1,x2,x3],np.float32)
  y3 = profile(x,lgm,ab,bc)/np.power(10.0,lgm)
  ax1.plot(r,y3, '-', c=(0,1.0,1.0,1.0),lw=2)

  x1 = r*0.0
  x2 = r*0.0
  x3 = r
  x = np.array([x1,x2,x3],np.float32)
  y3 = profile(x,lgm,ab,bc)/np.power(10.0,lgm)
  ax1.plot(r,y3, '-', c=(1.0,1.0,0.0,1.0),lw=2)




plt.show()
