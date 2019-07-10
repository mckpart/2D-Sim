import numpy as np
from matplotlib import pyplot as plt
import math
import yaml

def dist(x1,x2,y1,y2):
   return math.sqrt((x2-x1)**2 + (y2-y1)**2)

def density(): 			# compute the number of particles per area
   for n in range(step): 

      area = math.pi * (deltaR * (n+1))**2
      part = float(N[n])
      g[n] = part/area

def numberAtRadius():	# compute the number of particles =< a distance

   x1 = position[4,0]
   y1 = position[4,1]

   for n in range(step):

      num = 0

      for k in range(int(n_part_tot - 1)): 
       
         x2 = position[4,(k + 1) * 2]
         y2 = position[4,(k + 1) * 2 + 1]
         
         d = dist(x1,x2,y1,y2)
      
         if(d < deltaR * (n+1)): 
            num = num + 1

            print num,n
            N[n] = num	

######## read in .yaml parameters #######

with open("params.yaml",'r') as yf:
    yaml_dict = yaml.safe_load(yf)

n_part_tot  = yaml_dict["totalParticles"]

######### initialize lists and read in position data ######

deltaR = .03
position = []

file = open( "positions.txt", "r" )
for line in file:
   row = line.split()
   row = [float(i) for i in row]
   position.append(row)

position = np.asarray(position)
# print position

# compute stepsize, initialize number and density arrays #

step = int(2/deltaR)

N = [0] * step
N = np.asarray(N)

r = np.linspace(0 + deltaR,2, step)

g = [0] * step
g = [float(i) for i in g]
g = np.asarray(g)

numberAtRadius(); 
density(); 

# plot the the density from a specific particle as a function
# of radial distance from the particle

plt.plot(r,g)
plt.xlim([0,2])
plt.title("Radial Density Function")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.show()

