import matplotlib.pyplot as plt
import sys
import numpy as np

file = open(sys.argv[1])
v = np.array([])
for line in file:
    try:
        t = np.array([float(t) for t in line.split()])
        v = np.append(v, t)
    except: pass

v.resize(len(v)/3, 3)
v[:,0] += -np.amin(v[:,0])
plt.grid(True)
plt.vlines(v[:,0], [0], v[:,2], color='red', lw=2, label='+')
plt.vlines(v[:,0], [0], v[:,1], color='blue', lw=2, label='-')
#plt.legend(loc='upper right')
plt.xlim([-10, np.amax(v[:,0])+10])
#plt.xlim([-10, 100])
plt.show()

