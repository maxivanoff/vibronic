import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
majorLocator   = MultipleLocator(100)
majorFormatter = FormatStrFormatter('%d')
file = open(sys.argv[1])
v = np.array([])
data = {}
for line in file:
    for n in ['Vab', 'w+', 'w-', 'b+', 'b-']:
        if line.split()[0] == n: data[n] = float(line.split()[2])
    try:
        t = np.array([float(t) for t in line.split()])
        if len(t)==3:
            v = np.append(v, t)
    except: pass
v.resize(len(v)/3, 3)
v[:,0] += -np.amin(v[:,0])
fig, ax = plt.subplots()
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
plt.grid(True)
plt.title('Vab = %d \nw- = %d b- = %.2f\n w+ = %d b+ = %.2f' % (data['Vab'], data['w-'], data['b-'], data['w+'], data['b+']))
plt.vlines(v[:,0], [0], v[:,2], color='blue', lw=2, label='-')
plt.vlines(v[:,0], [0], v[:,1], color='red', lw=2, label='+')
plt.legend(loc='upper right')
plt.xlim([-10, np.amax(v[:,0])+10])
plt.xlim([-10, 1000])
plt.savefig("%.2f.pdf" % data['b-'])
plt.show()

