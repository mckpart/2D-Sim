import matplotlib.pyplot as plt
import numpy as np
import math

dens = np.linspace(.1,.9,9)
print dens
pressure = [.20545, .454075, .782773, 1.20705,1.83991, 2.73406, 4.58495 \
        , 7.67438, 13.3525]

for k in range(9):
    plt.plot(dens[k],pressure[k],'ko')

plt.ylabel('Pressure')
plt.xlabel('Density')

plt.xlim([0,max(dens) + .1])
plt.ylim([0,max(pressure) + 1])

plt.show()
