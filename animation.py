import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt 
from matplotlib.animation import FuncAnimation
import yaml
import math

def dist(x1,x2,y1,y2):
   return math.sqrt((x2-x1)**2 + (y2-y1)**2)

######## read in .yaml parameters #######

with open("params.yaml",'r') as yf:
   yaml_dict = yaml.safe_load(yf)

radius_1   = float(yaml_dict["particleRadius_1"])
radius_2   = float(yaml_dict["particleRadius_2"])     
n_part_1   = yaml_dict["type1_Particles"]
n_part_2   = yaml_dict["type2_Particles"]
n_part_tot = yaml_dict["totalParticles"]

# numIter    = yaml_dict["numberUpdates"]
rigidBC    = yaml_dict["rigidBoundary"]
c_linkers  = yaml_dict["crosslinkers"]

boxLength  = yaml_dict["boxLength"]
# periodBC    = yaml_dict["periodicBoundary"]

######### initialize lists and read in position data ##########

patches = []
position = []
x,y = [],[]

fig,ax = plt.subplots()

colors = mpl.cm.rainbow(np.linspace(0,1, n_part_tot))
for i in range(n_part_1):
    patches += [plt.Circle((0,0), radius_1, color= colors[i])]
for i in range(n_part_2):
    patches += [plt.Circle((0,0), radius_2, color = colors[i + n_part_1])]    

file = open( "positions.txt", "r" )
for line in file:
    row = line.split()
    row = [float(i) for i in row]
    position.append(row)

position = np.asarray(position)
# print position

numIter = len(position[:,0]) 
def init():

   if(c_linkers == 1 and rigidBC != 1):
      ax.set_xlim(-2,2)
      ax.set_ylim(-2,2)
   else:  
      ax.set_xlim(-1 * boxLength,boxLength)       
      ax.set_ylim(-1 * boxLength,boxLength)
        
   for i in range(n_part_tot):
      ax.add_patch(patches[i])

   return patches

def update(frame):

   frame = int(frame)
   print "frame is ", frame
   
   for k in range(n_part_tot):
      x = position[frame][k * 2]
      y = position[frame][(k * 2) + 1]

      patch = patches[k]
      patch.center = (x,y)
      patches[k] = patch

   return patches


t = np.linspace(0,numIter - 1,numIter) 

anim = FuncAnimation(fig, update, frames = t,init_func = init,blit = True)
plt.show()
