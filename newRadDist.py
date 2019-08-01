import math 
import matplotlib.pyplot as plt
import numpy as np
import yaml

def calculateRadDist(N, numberTrials,total_dens):
    G = np.zeros(len(N))
    G = map(lambda x: x/(n_particles * numberTrials * total_dens), N)

    return G

def calcDensity(num_dens, r_iter, deltaR):
    N = np.zeros(len(num_dens))
    for k in range(len(r_iter)): 
        if(k == 0):
            area = math.pi * deltaR**2 / 4.0
        else: 
            area = 2.0 * math.pi * r_iter[k]*deltaR
        N[k] = num_dens[k]/area    

    return N 
# Read in the number density
# The data has not yet been normalized 

with open('numDensity.txt','r') as f: 
    n_dens = f.read().split(' ')

#print n_dens
n_dens = [float(i) for i in n_dens if i != '']    
n_dens = np.asarray(n_dens)
#print n_dens

# read in system parameters from yaml file

with open("params.yaml",'r') as yf:
    yaml_dict = yaml.safe_load(yf)

boxL = yaml_dict["boxLength"]
n_particles = yaml_dict["totalParticles"]
sigma = yaml_dict["sigma"]
redDens = yaml_dict["reducedDens"]
n_eq = yaml_dict["equilibriate_sweep"]
n_updates = yaml_dict["numberUpdates"]
interval = yaml_dict["data_collect_interval"]

if(boxL == 0): 
    boxL = sigma * math.sqrt(n_particles/redDens)
if(sigma == 0):
    sigma = boxL * math.sqrt(redDens/n_particles)
if(redDens == 0):
    redDens = n_particles * sigma**2 / boxL**2

trunc_dist = 2.5 * sigma

total_dens = n_particles / boxL**2

print boxL

##### PREPARE SYSTEM ##############
numberTrials = (n_updates-n_eq)/interval    # this needs to be added to the yaml file
deltaR = sigma/20.0 

r_iter = np.linspace(0,.5*boxL,(.5 * boxL)/deltaR + 1.0)
#print r_iter
N = calcDensity(n_dens,r_iter,deltaR)
G = calculateRadDist(N,numberTrials,total_dens);
#print G

# plot the results with comparison values y = 1
# and x = 2^(1/6)

plt.plot(r_iter,G)
plt.axvline(x = 2.0**(1.0/6.0),color = 'r', linestyle = '--')
plt.axhline(y = 1, linestyle = '--')

plt.title('Radial Distribution')
plt.legend(['RDF',r'$2^{1/6}$'])
plt.xlabel(r'$\frac{r}{\sigma}$')
plt.ylabel(r'$g(\frac{r}{\sigma})$')

plt.xlim([0,boxL/2])
plt.show()
