import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import math
import yaml

def dist(x1,x2,y1,y2):
   return math.sqrt((x2-x1)**2 + (y2-y1)**2)

def calculate_g():	# compute the number of particles =< a distance

    sys_dens = n_part_tot / (math.pi * boxLength**2) # system density

    for n in range(step):  
      num = 0.0 
      r_curr = deltaR * n  # current distance

      area = math.pi * ((r_curr + deltaR)**2 - r_curr**2) # current area
       
      for count in range(n_positions): # averages over position and time
         for p in range(n_part_tot):

            x1 = position[count,2 * p]
            y1 = position[count,(2 * p) + 1]

            for k in range(int(n_part_tot)): 
               if(k != p):

                  x2 = position[count,k * 2]
                  y2 = position[count,(k * 2) + 1]
            
                  d = dist(x1,x2,y1,y2)
      
                  if(d > r_curr and d <= r_curr + deltaR): # counts the number
                     num = num + 1          # of particles at current distance
          
      avgNum = num /(.5 * n_positions * (n_part_tot )) # agrees with simple model
      g[n] = avgNum/(area * sys_dens)   # normalizes the function

######## read in .yaml parameters #######

with open("params.yaml",'r') as yf:
    yaml_dict = yaml.safe_load(yf)

n_part_tot  = yaml_dict["totalParticles"]
boxLength = yaml_dict["boxLength"]
restLength = yaml_dict["restLength"]
sigma = yaml_dict["sigma"]
radius_1 = yaml_dict["particleRadius_1"]
LJ = yaml_dict["lennardJones"]

######### initialize lists and read in position data ######

deltaR = boxLength/100.0
position = []

file = open( "positions.txt", "r" )
for line in file:
   row = line.split()
   row = [float(i) for i in row]
   position.append(row)

position = np.asarray(position)

n_positions = len(position[:,0])
print n_positions

# compute stepsize, initialize arrays #

step = int(boxLength/deltaR) + 1

r = np.linspace(0, boxLength, step)
g = np.zeros(step)

if(LJ == 1):        # reduced length for lennard jones system 
   r = [i/(sigma) for i in r]
   radius_1 = radius_1/sigma
   restLength = 2.0**(1.0/6.0)

#### generate densities and plot ###########3

calculate_g(); 
plt.plot(r,g)

plt.axvline(x = restLength, color = 'r', linestyle = '--')
plt.axvline(x = radius_1 * 2, color = 'g', linestyle = '--')

plt.xlim([0,max(r)])
plt.ylim([0,max(g) * 1.1])

plt.title("Radial Distribution Function")
plt.xlabel(r'$\frac{r}{\sigma}$')
plt.ylabel(r'$g(\frac{r}{\sigma})$')
plt.legend([r'$g(\frac{r}{\sigma})$', 'rest length', 'diameter'])

plt.show()
 
