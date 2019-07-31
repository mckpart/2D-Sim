import matplotlib.pyplot as plt
import numpy as np
import math

dens = [.7, .75, .80255,.82, .83, .84,.855,.9]
print dens
# pressure = [.212256, .463658, .777451, 1.25733, 1.92992, 2.94772,\
#          4.80433, 8.10383, 14.3627]

correction = [-.2356740339,-.2705441716,-.30978462,-0.32340249,-.331338453,-.3393706089, \
              -.3515992054,-.3895836072]
pressure_1 = [.982236, 1.45257, 2.21196,2.5772,2.79178,3.15011,3.4018,5.50457]
pressure_2 = [0.71,1.20,2.01,2.28,2.45,2.53,2.18,3.60]

for k in range(7):
    plt.plot(dens[k],pressure_1[k],'ko')
    plt.plot(dens[k],pressure_2[k], 'bX')
    plt.plot(dens[k],pressure_1[k] + correction[k],'ro')

plt.title(r'Isotherm $T^* = 0.7$')
plt.legend(['actual','expected','shifted'])
plt.ylabel('Pressure')
plt.xlabel('Density')

plt.xlim([.7,max(dens) + .1])
plt.ylim([0,max(pressure_1) + 1])

############### Phase diagram in (T,P) plane #####################

plt.figure(2)

T = [0.45, 0.55,0.7,.75,1]
P_act = [.135056, .468287, .98226, 1.13619,1.94988]
P_exp = [-.04,.23,.71,.92,1.67]
c = -.2356740339

for k in range(5):
    plt.plot(T[k],P_act[k],'ko')
    plt.plot(T[k],P_exp[k],'bX')
    plt.plot(T[k],P_act[k] + c,'ro')

plt.title(r'Phase plane of $(T^*,P^*)$ $\rho^* = 0.7$')
plt.legend(['actual','expected','shifted'])
plt.xlabel('Temperature')
plt.ylabel('Pressure')

plt.xlim([.4,1.1])
plt.ylim([-.1,2])
plt.show()
