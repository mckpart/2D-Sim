import math
import numpy as np
from matplotlib import pyplot as plt

def potentialEnergy(r):
    # return -1 * math.exp(-1 * (r - r0)**2) + 1
    return (r-r0)**2 * math.exp(-12 * (r - r0)**2)

r0 = .2
r = np.linspace(0,1,100) 
U = map(potentialEnergy,r)

plt.plot(r,U)
plt.xlim([0,1])
plt.show()
