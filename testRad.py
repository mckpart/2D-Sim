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

      # if(n * deltaR > boxLength):
      #    g[n] = 1
      # else:
         part = N[n]
         g[n] = part/(area * 0.5 * n_part_tot) 

def avgNumberAtRadius():	# compute the number of particles =< a distance

   # x1 = position[4,0]
   # y1 = position[4,1]

   for n in range(step): # THIS MAY BE TOO LARGE
      num = 0.0 
       
      for count in range(n_positions/2):
         c = 0
         for p in range(n_part_tot):

            x1 = position[2 * count,2 * p]
            y1 = position[2 * count,(2 * p) + 1]

            dist_x = boxLength - abs(x1)
            dist_y = boxLength - abs(y1)

            if(dist_x > dist_y):
                r_max = dist_x
            else:
                r_max = dist_y
            
            if(n * deltaR < r_max): 
               for k in range(int(n_part_tot)): 
                  if(k != p):
                     c = c + 1 
                     x2 = position[2 * count,k * 2]
                     y2 = position[2 * count,(k * 2) + 1]
            
                     d = dist(x1,x2,y1,y2)
      
                     if(d > deltaR * n and d <= deltaR * (n+1)): 
                        num = num + 1
            # else:
               # c = 1
               # num = 0
               # print "c is ", c,"rmax is ",r_max 
               # print "n * delR is ", n * deltaR
               # break
         if(c == 0): 
             avgNum = 0
         else:     
            avgNum = num / (2 * c  * n_positions)
         
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

deltaR = .02
position = []

file = open( "positions.txt", "r" )
for line in file:
   row = line.split()
   row = [float(i) for i in row]
   position.append(row)

position = np.asarray(position)
# print position

n_positions = len(position[:,0])
print n_positions

# compute stepsize, initialize number and density arrays #

graphLen = boxLength

step = int(graphLen/deltaR)
print "step is", step


sys_dens = n_part_tot / (2*boxLength)**2

N = [0] * step
N = [float(i) for i in N]
N = np.asarray(N)

r = np.linspace(0, graphLen, step)
r = [i/(2*radius_1) for i in r]
print "r is", r
g = [0] * step
g = [float(i) for i in g]
g = np.asarray(g)

avgNumberAtRadius(); 
density(); 

# G = interp1d(r,g, kind = 'cubic')
rnew = np.linspace(0,graphLen, step * 3)

plt.plot(r,g)
# plt.plot(rnew,G(rnew))

plt.axvline(x = restLength, color = 'r', linestyle = '--')
plt.axvline(x = radius_1 * 2, color = 'g', linestyle = '--')

plt.xlim([0,max(r)])
plt.ylim([0,max(g) * 1.1])

plt.title("Radial Density Function")
plt.xlabel("r/ $/sigma$")
plt.ylabel("g(r)")
plt.legend(['g(r)', 'rest length', 'diameter'])

plt.show()
 
