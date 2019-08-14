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

trunc = 2.5 * sigma

boxLength  = yaml_dict["boxLength"]

f = open('forces.txt','r')
vir = f.read().split(' ')
vir = [float(i) for i in vir if(i != '' and i != '\n')]

def pressure_corr_LJ(d): 
    t = trunc/sigma
    c = 6*math.pi* d**2 * (.8*(1/t)**10.0 - (1/t)**4.0)
    return c 

def calcPressure(d,t,n_part,vir):
   avg_virial = 0
   for k in range(len(vir)): 
       avg_virial = avg_virial + vir[k]
   avg_virial = avg_virial/len(vir)
   p = d * (t + avg_virial/(2 * n_part))
   if(interact_type == 1): 
       corr = pressure_corr_LJ(d)
   else:
       corr = 0
   return p + corr 

pressure = calcPressure(red_dens,red_temp,n_part,vir)
print("reduced pressure: ",pressure)
