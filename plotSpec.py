import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
majorLocator   = MultipleLocator(200)
majorFormatter = FormatStrFormatter('%d')
file = open(sys.argv[1])
v = np.array([])
data = {}
for line in file:
    for n in ['Vab', 'w+', 'w-', 'b+', 'b-', 'w', 'b']:
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
d = v[:,0][np.argmax(v[:,1])] - v[:,0][np.argmax(v[:,2])]
try:
    plt.title(r'Vab = %.3f $\delta$ = %.4f\nw- = %d b- = %.2f\n w+ = %d b+ = %.2f' % (data['Vab'],d, data['w-'], data['b-'], data['w+'], data['b+']))
except: pass
try:
    plt.title(r'Vab = %.3f $\delta$ = %.4f' % (data['Vab'], d) + '\n w = %d b = %.2f' % ( data['w'], data['b']))
except: pass
plt.vlines(v[:,0], [0], v[:,2], color='blue', lw=2, label='-')
plt.vlines(v[:,0], [0], v[:,1], color='red', lw=2, label='+')
n = len(v[:,1])
plt.legend(loc='upper right')
plt.xlim([-10, np.amax(v[:,0])+10])
plt.xlim([-10, 1200])
plt.savefig("%.2f.pdf" % data['b'])
plt.show()

