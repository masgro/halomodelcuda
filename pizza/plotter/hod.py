import matplotlib.pyplot as plt
import pylab as pp
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import rcParams
from scipy.special import erf

NCEN_LOGMMIN = 12.481
NCEN_SIGMA   =  0.15
NSAT_LOGM0   = 13.129
NSAT_LOGM1   = 13.661
NSAT_ALPHA   =  1.0

def ncen(m):
  global NCEN_LOGMMIN, NCEN_SIGMA
  y = 0.5*(1.0 + erf((m - NCEN_LOGMMIN)/NCEN_SIGMA))
  return y

def nsat(m):
  global NSAT_LOGM0, NSAT_LOGM1, NSAT_ALPHA
  y = np.power(10.,m) - np.power(10.,NSAT_LOGM0)
  y = y/np.power(10.,NSAT_LOGM1)
  y = np.power(y,NSAT_ALPHA)
  return y

def hod(m):
  y = ncen(m) + nsat(m)
  return y

f = plt.figure(figsize=(10,10))
ax = f.add_subplot(111)

x = np.linspace(9,16)
y = hod(x)
ax.plot(x,np.log10(y))

NCEN_LOGMMIN = 12.281
NCEN_SIGMA   =  0.15
NSAT_LOGM0   = 13.129
NSAT_LOGM1   = 13.861
NSAT_ALPHA   =  1.0
x = np.linspace(9,16)
y = hod(x)
ax.plot(x,np.log10(y))

plt.show()
