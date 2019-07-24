import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as tck

f_file = open("forces.txt","r")
forces = f_file.read().split(' ')
forces = [float(i) for i in forces if i != '']

e_file = open("energies.txt","r")
energy = e_file.read().split(' '); print energy; 
energy = [float(i) for i in energy if i != '']

print forces

numIter = np.linspace(1,len(forces),len(forces))
print numIter

plt.subplot(211)
plt.plot(numIter,forces)
plt.xlabel('Number of Sweeps')
plt.ylabel(r'$F \cdot r$')

ax = plt.gca()
ticks = tck.FuncFormatter(lambda x, pos: '{0:g}'.format(x*10))
ax.xaxis.set_major_formatter(ticks)

plt.subplot(212)
plt.plot(numIter,energy)
plt.xlabel('Number of Sweeps')
plt.ylabel('Free energy')

ax = plt.gca()
ticks = tck.FuncFormatter(lambda x, pos: '{0:g}'.format(x*10))
ax.xaxis.set_major_formatter(ticks)

plt.show()
