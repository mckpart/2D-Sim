import matplotlib.pyplot as plt
import numpy as np
import math

dens = np.linspace(.1,.9,9)
print dens
pressure = [.212256, .463658, .777451, 1.25733, 1.92992, 2.94772,\
         4.80433, 8.10383, 14.3627]

for k in range(9):
    plt.plot(dens[k],pressure[k],'ko')

plt.title('Pressure v. Density')
plt.ylabel('Pressure')
plt.xlabel('Density')

plt.xlim([0,max(dens) + .1])
plt.ylim([0,max(pressure) + 1])

plt.show()
