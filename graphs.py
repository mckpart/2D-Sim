import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as tck

# f_file = open("forces.txt","r")
# forces = f_file.read().split(' ')
# forces = [float(i) for i in forces if i != '']
# 
# e_file = open("energies.txt","r")
# energy = e_file.read().split(' '); print energy; 
# energy = [float(i) for i in energy if i != '']
# 
# print forces
# 
# numIter = np.linspace(1,len(forces),len(forces))
# print numIter
# 
# plt.subplot(211)
# 
# plt.plot(numIter,forces)
# 
# plt.title("Square Lattice")
# plt.xlabel('Number of Sweeps')
# plt.ylabel(r'$F \cdot r$')
# 
# ax = plt.gca()
# ticks = tck.FuncFormatter(lambda x, pos: '{0:g}'.format(x*20))
# ax.xaxis.set_major_formatter(ticks)
# 
# plt.subplot(212)
# plt.plot(numIter,energy)
# plt.xlabel('Number of Sweeps')
# plt.ylabel('Free energy')
# 
# ax = plt.gca()
# ticks = tck.FuncFormatter(lambda x, pos: '{0:g}'.format(x*20))
# ax.xaxis.set_major_formatter(ticks)
# 
# plt.show()
a = 1.0
k = 1.0
r0 = 0.0

x = np.linspace(-10,10,10000)
U = [a*k/2 *(r-r0)**2 * math.exp(-k/2 *(r-r0)**2) for r in x]
der_func = lambda r: -a*k*(r-r0)*math.exp(-k/2*(r-r0)**2)*(1-k/2*(r-r0)**2)
print(der_func(0))
d = [der_func(i) for i in x]

plt.plot(x,U)
plt.plot(x,d)
plt.axhline(y = 0) 
plt.show()
