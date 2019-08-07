import math 
import matplotlib.pyplot as plt
import numpy as np
import yaml

# Normalizes the set based on the total 
#   number of trials in the set
# Returns the radial distribution function
def calculateRadDist(N,numberTrials,total_dens):
    G = np.zeros(len(N))
    G = [x/(n_particles*numberTrials*total_dens) for x in N]
    print(n_particles,len(N),total_dens)
    return G

# Creates an unnormalized vector of densities
# Accepts a vector of the number density and
#   the value of deltaR
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
def readData(file_name): 
    with open(file_name,'r') as f: 
        n = f.read().split(' ')

        n = [float(i) for i in n if i != '']    
        n = np.asarray(n)
    f.close()
    return n 

n_dens = readData('numDensity.txt')
par_dens = readData('par_numDensity.txt')
antp_dens = readData('antp_numDensity.txt')

dens_vec = [n_dens,par_dens,antp_dens,par_dens + antp_dens]
print(dens_vec)
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
red_temp = yaml_dict["reducedTemp"]

if(boxL == 0): 
    boxL = sigma * math.sqrt(n_particles/redDens)
if(sigma == 0):
    sigma = boxL * math.sqrt(redDens/n_particles)
if(redDens == 0):
    redDens = n_particles * sigma**2 / boxL**2

trunc_dist = 2.5 * sigma
total_dens = n_particles / boxL**2

print(boxL)

##### PREPARE SYSTEM ##############
numberTrials = (n_updates-n_eq)/interval - 1  # this needs to be added to the yaml file
deltaR = sigma/20.0 
r_iter = np.linspace(0,.5*boxL,(.5 * boxL)/deltaR + 1.0)

##### PLOT DATA ###################
choice = input("Do you want all of the plots?(y/n)")

# 'y' will create an image with the regular RDF, parallel RDF, and antiparallel RDF
# anything else will create an image with the regular RDF
if(choice == 'y'):  
    fig, axs = plt.subplots(1,3, figsize=(12, 3), 
               facecolor='w', edgecolor='k',squeeze = False)
    fig.subplots_adjust(hspace = .5, wspace=.2)
    num = 3
else: 
    fig, axs = plt.subplots(figsize = (7,5),
               facecolor = 'w',edgecolor='k',squeeze = False)
    num = 1

titles = ['Radial Distribution Function','Parallel RDF','Antiparallel RDF','a + p']
for k in range(num):
    N = calcDensity(dens_vec[k],r_iter,deltaR)
    G = calculateRadDist(N,numberTrials,total_dens)
    
    axs[0,k].plot(r_iter,G)
    axs[0,k].axvline(x = 2.0**(1.0/6.0),color = 'r', linestyle = '--')
    axs[0,k].axhline(y = 1, linestyle = '--')

    axs[0,k].set_title(titles[k])
    axs[0,k].legend(['RDF',r'$2^{1/6}$'])
    axs[0,k].set_xlabel(r'$\frac{r}{\sigma}$')

    axs[0,k].set_xlim([0,boxL/2 - 2 * deltaR])
    if(k == 0):
        axs[0,k].set_ylabel(r'$g(\frac{r}{\sigma})$')
        max_y = max(G) * 1.1
    axs[0,k].set_ylim([0,max_y])


### GIVES THE OPTION TO SAVE THE GENERATED IMAGE ##
s = input("would you like to save the image? (y/n)")
if(s == 'y'):
    file_name = 'RDF_dens_' + str(redDens).replace('.','_') + \
                'temp_' + str(red_temp).replace('.','_') + \
                '.png'
    
    print(file_name)
    plt.savefig('data/' + file_name)
    plt.close()

plt.show()   
