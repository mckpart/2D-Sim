import math
import yaml

with open("params.yaml",'r') as yf:
   yaml_dict = yaml.safe_load(yf)

radius   = yaml_dict["particleRadius"]
n_part = yaml_dict["totalParticles"]

sigma = yaml_dict["sigma"]
interact_type = yaml_dict["interactionType"]
red_dens = yaml_dict["reducedDens"]
red_temp = yaml_dict["reducedTemp"]

# rigidBC    = yaml_dict["rigidBoundary"]
# periodBC    = yaml_dict["periodicBoundary"]

boxLength  = yaml_dict["boxLength"]

f = open('forces.txt','r')
vir = f.read().split(' ')
vir = [float(i) for i in vir if(i != '' and i != '\n')]

def calcPressure(d,t,n_part,vir):
   avg_virial = 0
   for k in range(len(vir)): 
       avg_virial = avg_virial + vir[k]
   avg_virial = avg_virial/len(vir)
   p = d * (t + avg_virial/(2 * n_part))
   return p

pressure = calcPressure(red_dens,red_temp,n_part,vir)
print("reduced pressure: ",pressure)
