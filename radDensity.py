import numpy as np
from scipy.ndimage.filters import gaussian_filter1d
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt
import math
import yaml

def dist(x1,x2,y1,y2):
   return math.sqrt((x2-x1)**2 + (y2-y1)**2)

def density(): 			# compute the number of particles per area
   for n in range(step): 

      if(n == 0): 
         area = math.pi * deltaR**2
      else:     
         area = 2 * math.pi * (n * deltaR) * deltaR

      part = N[n]
      g[n] = part/(area * n_part_tot/(2 * boxLength)**2)

def avgNumberAtRadius():	# compute the number of particles =< a distance

   # x1 = position[4,0]
   # y1 = position[4,1]

   for n in range(step):
      num = 0.0

      for p in range(n_part_tot):

         x1 = position[10,2 * p]
         y1 = position[10,(2 * p) + 1]

         for k in range(int(n_part_tot)): 
            if(k != p): 
               x2 = position[10,k * 2]
               y2 = position[10,(k * 2) + 1]
         
               d = dist(x1,x2,y1,y2)
      
               if(d > deltaR * n and d <= deltaR * (n+1)): 
                  num = num + 1	

      avgNum = num / (2 * n_part_tot)
      N[n] = avgNum 
      # print "the average number is:", avgNum 

######## read in .yaml parameters #######

with open("params.yaml",'r') as yf:
    yaml_dict = yaml.safe_load(yf)

n_part_tot  = yaml_dict["totalParticles"]
boxLength = yaml_dict["boxLength"]
restLength = yaml_dict["restLength"]
radius_1 = yaml_dict["particleRadius_1"]

######### initialize lists and read in position data ######

deltaR = .01
position = []

file = open( "positions.txt", "r" )
for line in file:
   row = line.split()
   row = [float(i) for i in row]
   position.append(row)

position = np.asarray(position)
# print position

# compute stepsize, initialize number and density arrays #

step = int(boxLength/deltaR)

N = [0] * step
N = [float(i) for i in N]
N = np.asarray(N)

r = np.linspace(0,boxLength, step)

g = [0] * step
g = [float(i) for i in g]
g = np.asarray(g)

avgNumberAtRadius(); 
density(); 

G = interp1d(r,g, kind = 'cubic')
rnew = np.linspace(0,boxLength, step * 3)

plt.plot(r,g)
plt.plot(rnew,G(rnew))

plt.axvline(x = restLength, color = 'r', linestyle = '--')
plt.axvline(x = radius_1 * 2, color = 'g', linestyle = '--')

plt.xlim([0,boxLength])
plt.ylim([0,max(G(rnew)) * 1.1])

plt.title("Radial Density Function")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.legend(['g(r)', 'rest length', 'radius'])

plt.show()
 
