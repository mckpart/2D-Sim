from shutil import copyfile
from shutil import move
#from sys import exit
import os
import yaml

name = input("please enter folder name: ")

pos_curr = "positions.txt"
rad_curr = "radialDistance.txt"
eng_curr = "energies.txt"
vir_curr = "forces.txt"
yaml_curr = "params.yaml"

old = [pos_curr, rad_curr, eng_curr, vir_curr, yaml_curr]

dir_name = "data/" + name
if not os.path.exists(dir_name):
    os.makedirs(dir_name)
else:
    print("The file ",name," already exists.")
pos_copy = "pos_" + name + ".txt"
rad_copy = "r_dist_" + name + ".txt"
eng_copy = "energy_" + name + ".txt"
vir_copy = "vir_" + name + ".txt"
yaml_copy = name + ".yaml"

new = [pos_copy,rad_copy,eng_copy,vir_copy,yaml_copy]
#print(new)

for k in range(5):
    copyfile(old[k],new[k])
    move(new[k], dir_name)
