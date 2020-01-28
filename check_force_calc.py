import math
import numpy as np

ff = open('forces.txt')
vir = ff.read().split(' ')
vir = [float(i)/30.0 for i in vir if i != '']
vir = np.asarray(vir)

f = open('avgForcePerParticle.txt','r')
run = f.read().split('\n')

M = []
for i in range(len(run)):
    temp = run[i].split(' ')
    temp = [float(i) for i in temp if i != '']
    temp = np.asarray(temp)
    M.append(temp)
#     print temp
print "M = ", M 
sum_vir = np.zeros(len(M))
for i in range(len(M)): 
    v = M[i]
    sum_curr = 0
    for k in range(len(v)/2): 
        x = v[2*k]
        y = v[2*k+1]

        sum_curr = sum_curr + x + y
#     print(len(v))
    sum_vir[i] = sum_curr * 1.0 / 30.0 
other_vir = 0
for i in range(len(vir)): 
    other_vir = other_vir + vir[i]
avg_vir = 0
for i in range(len(sum_vir)): 
    avg_vir = avg_vir + sum_vir[i]

other_vir = other_vir / len(vir)
avg_vir = avg_vir / len(sum_vir)
print avg_vir,other_vir

# print M
