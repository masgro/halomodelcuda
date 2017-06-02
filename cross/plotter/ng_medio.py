import matplotlib.pyplot as plt
import numpy as np

filename = "../integrando_ng_medio.dat"
x,y,z = np.genfromtxt(filename,unpack=True,usecols=(0,1,2))

plt.plot(x,y)

##data=np.genfromtxt("ho6.dat",unpack=True)
##ii=np.where(data[2]>0)
##y_meas=data[3][ii]/data[2][ii]
##x=np.log10(data[0][ii])
##plt.plot(x,y_meas,'rx',markersize=12)
plt.show()
