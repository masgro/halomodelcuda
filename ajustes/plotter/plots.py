import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

filename = 'Nu_M.dat'
nu,mo,mi = np.genfromtxt(filename,unpack=True)

gs = gridspec.GridSpec(3,1)

f = plt.figure()
ax1 = f.add_subplot(gs[0:2,0])
ax2 = f.add_subplot(gs[2,0],sharex=ax1)

ax1.plot(nu,mo)
ax1.plot(nu,mi)

ax2.plot(nu,mi-mo)
ax2.set_ylim(-0.01,0.01)
ax1.set_xlim(9.0,16.0)
ax1.set_ylim(-1.0,1.0)

plt.show()

