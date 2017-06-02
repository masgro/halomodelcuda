import matplotlib.pyplot as plt
import numpy as np

phi = 0.0
theta = 0.0

def set_psi(ab,bc,theta,phi):

  t = ab*ab
  t = 1./t
  s = bc*bc
  s = 1./s

  ct = np.cos(theta)
  ct2 = ct*ct
  st2 = 1. - ct2
  st = np.sqrt(st2)

  cf = np.cos(phi)
  sf = np.sin(phi)
  sf2 = sf*sf
  cf2 = cf*cf

  A = ct2*(sf2*t*s + cf2*s) + st2*t*s*s
  B = ct*np.sin(2.*phi)*(t*s - s)
  C = (sf2*s + cf2*t*s)

  psi = 0.5*np.arctan2(B,(A-C))

  if(((A - C)*np.cos(2.*psi) + B*np.sin(2.*psi)) > 0.0):
    psi = psi - np.pi*0.5

  return(psi)

def rotateA(X,phi,theta):
  ##A = [[-np.sin(phi), -np.cos(phi)*np.cos(theta), np.cos(phi)*np.sin(theta)],
  ##     [ np.cos(phi), -np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta)],
  ##     [           0,              np.sin(theta), np.cos(theta)]]


  A = [[ np.cos(phi), -np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta)], 
       [ np.sin(phi), np.cos(phi)*np.cos(theta), -np.cos(phi)*np.sin(theta)], 
       [0, np.sin(theta), np.cos(theta)]]

  Y = np.dot(A,X)
  return(Y)

def rotateB(X,phi,theta,psi):
  B = [[-np.sin(phi)*np.cos(psi)-np.cos(phi)*np.cos(theta)*np.sin(psi), np.sin(phi)*np.sin(psi)-np.cos(phi)*np.cos(theta)*np.cos(psi), np.cos(phi)*np.sin(theta)],
      [np.cos(phi)*np.cos(psi)-np.sin(phi)*np.cos(theta)*np.sin(psi), -np.cos(phi)*np.sin(psi)-np.sin(phi)*np.cos(theta)*np.cos(psi), np.sin(phi)*np.sin(theta)],
      [np.sin(theta)*np.sin(psi), np.sin(theta)*np.cos(psi), np.cos(theta)]]
  Y = np.dot(B,X)
  return(Y)



ab = 0.5
bc = 0.5

theta = np.random.rand()
theta = np.arccos(theta)
phi = np.random.rand()*2.0*np.pi

theta = np.pi*0.5
phi = 0.0

psi = set_psi(ab,bc,theta,phi)


x = 1.0
y = 1.0
z = 0.1

X = np.asarray((x,y,z))
Y = rotateA(X,phi,theta)

print X
print Y



#*f = st*st*(cf2*t*s+s*sf2) + ct2;
#
#r[0] = -sf*xt-cf*ct*yt+cf*st*zt; 
#r[1] = cf*xt-sf*ct*yt+sf*st*zt;
#r[2] = st*yt+ct*zt;

