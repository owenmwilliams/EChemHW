from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math

#Raw input setup
l = int(raw_input('Number of time steps.'))
DMA = float(raw_input('Dimensionless diffusion coefficient.'))

#Array setup          
CAold =[]
CAnew=[]
conc_final = []
current_cott = []
current = []
time = []

#Bulk concentration
for x in range(l+1):
    CAold.append(1.0)


#Iterative Calculations for Current/Concentration
k = 1
while k < l+1:
    conc_final.append(CAold)
    Z = CAold[0] + DMA*(CAold[1]-CAold[0])
    Z_sim = math.sqrt(l/DMA)*Z
    current.append(Z_sim)
    t = float(k)/float(l)
    time.append(t)
    
    T = float(k-0.5)/float(l)
    Zcott = 1/math.sqrt(3.141592*T)
    current_cott.append(Zcott)
    
    CAnew.append(0)
    for i in range(1,l):
        CAnew.append(CAold[i]+DMA*(CAold[i+1]-2*CAold[i]+CAold[i-1]))
    CAnew.append(1)
    CAold = CAnew
    CAnew = []
    k = k+1

#Concentration profile plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(len(conc_final)):
    temp_dist = []
    temp_time = []
    for j in range(len(conc_final[i])):
        temp_dist.append(j)
        temp_time.append(i)
    ax.scatter(temp_dist, temp_time, conc_final[i])
ax.set_xlabel('$Distance (arb.)$')
ax.set_ylabel('$Time (arb.)$')
ax.set_zlabel('$Concentration (arb.)$')
plt.title('$Concentration$ $Profile$ $;$ $t_k=$ $' + str(l) + '$ $;$ $DMA=$ $' + str(DMA) + '$')
plt.show()

#Current difference plot
diff = []
for i in range(len(current)):
    diff.append(current[i]/current_cott[i])
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time,diff,c='g',marker='s')
ax.set_xlabel('$t/t_k$', fontsize="large")
ax.set_ylabel('$Z_{sim}/Z_{Cottrell}$')
plt.title('$Difference$ $Between$ $Simulated$ $Current$ $and$ $Cottrell$ $Current$')
plt.show()

#Current plots
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(time,current_cott,c='b',marker='o',label='Cottrell Current')
ax.scatter(time,current,c='r',marker='^', label='Simulated Current')
ax.set_xlabel('$t/t_k$', fontsize="large")
ax.set_ylabel('$Current(arb.)$')
plt.title('$Simulated$ $Current$ $and$ $Cottrell$ $Current$')
ax.legend()
plt.show()
