import numpy as np
from matplotlib import cm
import matplotlib as mpl
from matplotlib import pyplot as plt 
import matplotlib.animation as ani_obj
from matplotlib.animation import FuncAnimation
import yaml
import math

def dist(x1,x2,y1,y2):
   return math.sqrt((x2-x1)**2 + (y2-y1)**2)

######## read in .yaml parameters #######

with open("params.yaml",'r') as yf:
   yaml_dict = yaml.safe_load(yf)

radius   = float(yaml_dict["particleRadius"])
# radius_2   = float(yaml_dict["particleRadius_2"])     
n_part_1   = yaml_dict["type1_Particles"]
n_part_2   = yaml_dict["type2_Particles"]
n_part_tot = yaml_dict["totalParticles"]

sigma = float(yaml_dict["sigma"])
interact_type = yaml_dict["interactionType"]
bound_type = yaml_dict["boundaryType"]
redDens = yaml_dict["reducedDens"]

# numIter    = yaml_dict["numberUpdates"]
# rigidBC    = yaml_dict["rigidBoundary"]
# periodBC    = yaml_dict["periodicBoundary"]
# c_linkers  = yaml_dict["crosslinkers"]

boxLength  = yaml_dict["boxLength"]
pos_file   = yaml_dict["animationFile"]

######### initialize lists and read in position data ##########

if(interact_type != 0):
    if(sigma == 0):
        sigma = boxLength * math.sqrt(redDens/n_part_tot)
    elif(boxLength == 0):
        boxLength = sigma * math.sqrt(n_part_tot/redDens)
    radius = 0.5 * sigma
patches = []
position = []
x,y = [],[]

fig,ax = plt.subplots()

file_1 = open('particle_type.txt','r')
types = file_1.read().split(' ')
types = [float(i) for i in types if i != '']
print(types)

for k in range(n_part_tot):
    if(types[k] == 1):
        c = 'red'
    else:
        c = 'blue'
    patches += [plt.Circle((0,0), radius, color = c)]
#color_1 = plt.cm.winter(np.linspace(0,1,n_part_1))
#color_2 = plt.cm.autumn(np.linspace(0,1,n_part_2))
#for i in range(n_part_1):
#    patches += [plt.Circle((0,0), radius, color = color_1[i])]
#for i in range(n_part_2):
#    patches += [plt.Circle((0,0), radius, color = color_2[i])]    

file = open(pos_file, "r" )
for line in file:
    row = line.split()
    row = [float(i) for i in row]
    position.append(row)

position = np.asarray(position)
print(position)

numIter = len(position[:,0]) 
print(numIter)

Writer = ani_obj.writers['ffmpeg']
writer = ani_obj.FFMpegWriter(fps=4, metadata=dict(artist='Me'), bitrate=1800)

def init():

   if(bound_type == 2):
      ax.set_xlim(-1 * boxLength,boxLength)
      ax.set_ylim(-1 * boxLength,boxLength)
   else:  
      ax.set_xlim(-.5 * boxLength,.5 * boxLength)       
      ax.set_ylim(-.5 * boxLength,.5 * boxLength)
        
   for i in range(n_part_tot):
      ax.add_patch(patches[i])

   return patches

def update(frame):

   frame = int(frame)
   print("frame is ", frame)
   
   for k in range(n_part_tot):
      x = position[frame][k * 2]
      y = position[frame][(k * 2) + 1]

      patch = patches[k]
      patch.center = (x,y)
      patches[k] = patch

   return patches


t = np.linspace(0,numIter - 1,numIter) 

anim = FuncAnimation(fig, update, frames = t,init_func = init,blit = True)

i = input("Would you like to save the movie?")
if(i == 'y'):
    anim.save('trial.mp4',writer=writer)
else:
    print("movie not saved")
plt.show()

