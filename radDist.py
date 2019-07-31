import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import math
import yaml

def dist(x1,x2):
   return abs(x2 - x1)

def calculate_g():	# compute the number of particles =< a distance

    sys_dens = n_part_tot / ((2 * boxLength)**2) # system density
    truncDist = 2.5 * sigma
    fact1 = 1
    fact2 = 1
    if(n_positions > 900): 
        fact1 = 50

    if(n_part_tot > 400 or n_positions > 400): 
        fact2 = 2
    elif(n_part_tot >= 100): 
        fact = 3
    elif(n_part_tot > 60): 
        fact = 2

    for n in range(step):  
      num = 0.0 
      r_curr = deltaR * n  # current distance

      area = math.pi * ((r_curr + deltaR)**2 - r_curr**2) # current area
       
      for count in range(n_positions/fact1): # averages over position and time
         for p in range(n_part_tot/fact2):

            x_curr = position[fact1 * count,2 * p * fact2]
            y_curr = position[fact1 * count,(2 * p * fact2) + 1]

            d_wall_curr_x = boxLength - abs(x_curr)
            d_wall_curr_y = boxLength - abs(y_curr)

            for k in range(int(n_part_tot/fact2)): 
               if(k != p):

                  x_comp = position[fact1 * count,k * 2 * fact2]
                  y_comp = position[fact1 * count,(k * 2 * fact2) + 1]
                  
                  x_dist = dist(x_curr,x_comp)
                  y_dist = dist(y_curr,y_comp)

                  d_wall_comp_x = boxLength - abs(x_comp)
                  d_wall_comp_y = boxLength - abs(y_comp)

                  if(d_wall_curr_x + d_wall_comp_x < truncDist and x_comp * x_curr < 0): 
                     x_dist = 2 * boxLength - x_dist

                  if(d_wall_curr_y + d_wall_comp_y < truncDist and y_comp * y_curr < 0):    
                     y_dist = 2 * boxLength - y_dist
                  
                  dist_tot = math.sqrt(x_dist**2 + y_dist**2)
      
                  if(dist_tot > r_curr and dist_tot <= r_curr + deltaR): # counts the number
                     num = num + 1          # of particles at current distance
          
      avgNum = float(fact1) * float(fact2)**2 * num /(0.5 * n_positions * (n_part_tot - 1 ))# agrees with simple model
      g[n] = avgNum/(area * sys_dens)   # normalizes the function

######## read in .yaml parameters #######

with open("params.yaml",'r') as yf:
    yaml_dict = yaml.safe_load(yf)

n_updates = yaml_dict["numberUpdates"]
n_part_tot  = yaml_dict["totalParticles"]
# beta = yaml_dict["beta"]
boxLength = yaml_dict["boxLength"]
# restLength = yaml_dict["restLength"]
sigma = yaml_dict["sigma"]
radius_1 = yaml_dict["particleRadius_1"]
LJ = yaml_dict["lennardJones"]
redDens = yaml_dict["reducedDens"]

######### initialize lists and read in position data ######

if(boxLength == 0):
    boxLength = sigma * math.sqrt(n_part_tot/redDens)
print "the box length is", boxLength
deltaR = float(0.5 * boxLength)/float(n_part_tot/2.0)
position = []
print("deltaR is ",deltaR)
file = open( "positions.txt", "r" )
for line in file:
   row = line.split()
   row = [float(i) for i in row]
   position.append(row)

position = np.asarray(position)

n_positions = len(position[:,0])
print n_positions

# compute stepsize, initialize arrays #

step = int(0.5 * boxLength/deltaR) + 1

r = np.linspace(0, 0.5 * boxLength, step)
g = np.zeros(step)

if(LJ == 1):        # reduced length for lennard jones system 
   r = [i/(sigma) for i in r]
#   radius_1 = radius_1/sigma
   restLength = 2.0**(1.0/6.0)

#### generate densities and plot ###########3

calculate_g();

plt.figure(figsize = (7.5,6.5))
plt.plot(r,g)

plt.axvline(x = restLength, color = 'r', linestyle = '--')
plt.axvline(x = radius_1 * 2, color = 'g', linestyle = '--')
plt.axhline(y = 1, linestyle = '--')

plt.xlim([0,max(r)])
plt.ylim([0,max(g) * 1.1])

plt.title("Radial Distribution Function")
plt.xlabel(r'$\frac{r}{\sigma}$')
plt.ylabel(r'$g(\frac{r}{\sigma})$')
plt.legend([r'$g(\frac{r}{\sigma})$', 'rest length', 'diameter'])

txt = "Iterations: " + str(n_updates) + "  Particles: " + str(n_part_tot) + \
        "  Box Length: " + str(boxLength) \
        + r'  $\sigma: $' + str(sigma) + r'  $\rho^*: $' + str(redDens)
plt.figtext(.5,.013,txt,wrap = True, ha = 'center')

plt.show()
