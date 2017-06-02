import numpy as np
import matplotlib.pyplot as plt

filename = '../integrando_ng_medio.dat'

data = np.genfromtxt(filename,unpack=True)

indx = (data[0] >= 10.82)
print np.trapz(data[1][indx],data[0][indx])

plt.plot(data[0],data[1])
plt.show()
