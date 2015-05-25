from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math

#Raw input setup
l = int(raw_input('Number of time steps.'))
DMA = float(raw_input('Dimensionless diffusion coefficient of species A.'))
DMB = float(raw_input('Dimensionless diffusion coefficient of species B.'))
IniA = float(raw_input('Molar concentration of species A.'))
IniB = float(raw_input('Molar concentration of species B.'))
k_o = float(raw_input('Electrode rate constant.'))
alpha = float(raw_input('Diffusion constant, alpha.'))
E_norm = float(raw_input('Potential step.'))
D_a = float(raw_input('Diffusion coefficient of A.'))
t_k = float(raw_input('Time constant.'))

#Array setup          
CAold =[]
CAnew=[]
CBold = []
CBnew =[]
conc_final_A = []
conc_final_B = []
current_cott = []
current = []
time = []

#Bulk concentration
for x in range(l+1):
    CAold.append(IniA)
    CBold.append(IniB)


#Iterative Calculations for Current/Concentration
k = 1
while k < l+1:
    conc_final_A.append(CAold)
    conc_final_B.append(CBold)

    t = float(k)/float(l)*t_k
    time.append(t)
    
    
    #Calculate current
    Z = k_o*math.sqrt((t_k)/(D_a)*math.exp(-alpha*E_norm*38.92)*CAold[0] - (k_o*math.sqrt((t_k)/(D_a)))*math.exp((1-alpha)*E_norm*38.92)*CBold[0])
    current.append(Z*DMA*96485.3)

    #Electrode reaction
    CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+Z*math.sqrt(DMA/l))
    CAnew.append(CAold[0]+DMA*(CAold[1]-CAold[0])-Z*math.sqrt(DMA/l))
    
    #Diffusion beyond the first box
    for i in range(1,l):
        CAnew.append(CAold[i]+DMA*(CAold[i+1]-2*CAold[i]+CAold[i-1]))
        CBnew.append(CBold[i]+DMB*(CBold[i+1]-2*CBold[i]+CBold[i-1]))

    CAnew.append(IniA)
    CBnew.append(IniB)
    CAold = CAnew
    CAnew = []
    CBold = CBnew
    CBnew = []
    k = k+1

#Concentration profile plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(len(conc_final_A)):
    temp_dist = []
    temp_time = []
    for j in range(len(conc_final_A[i])):
        temp_dist.append(j)
        temp_time.append(i)
    ax.scatter(temp_dist, temp_time, conc_final_A[i],c='b')
for i in range(len(conc_final_B)):
    temp_dist = []
    temp_time = []
    for j in range(len(conc_final_B[i])):
        temp_dist.append(j)
        temp_time.append(i)
    ax.scatter(temp_dist, temp_time, conc_final_B[i],c='r')
ax.set_xlabel('$Distance (arb.)$')
ax.set_ylabel('$Time (arb.)$')
ax.set_zlabel('$Concentration (M)$')
plt.title('$Concentration$ $Profile$ $;$ $t_k=$ $$' + str(l) + '$$ $;$ $(E-E^0)=$' + str(E_norm) + '$$ $;$ $k^0=$' + str(k_o))
plt.show()

#Current plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(time,current,c='r',marker='^', label='Simulated Current')
ax.set_xlabel('$t/t_k$', fontsize="large")
ax.set_ylabel('$Current(arb.)$')
plt.title('$Simulated$ $Current$')
ax.legend()
plt.show()
