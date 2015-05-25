from __future__ import division
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math


#Raw input setup
l = float(2500)
DMA = float(0.4)
DMB = float(0.4)
Da = float(0.00001)
IniA = float(1)
IniB = float(0)
k_o = float(0.001)
alpha = float(0.5)
E_o = float(0)
k_rxn_list = [1,0.1,0.01,0.001,0.0001]

master_time = []
master_current = []
for i in range(len(k_rxn_list)):
    k_rxn = k_rxn_list[i]
    scan = 0.1
    E_start = 0.15
    E_rev = -0.15
    #Array setup          
    CAold =[]
    CAnew=[]
    CBold = []
    CBnew =[]
    current = []
    time = []

    #Bulk concentration
    for x in range(int((l)+1)):
        CAold.append(IniA)
        CBold.append(IniB)


    #Iterative Calculations for Current/Concentration

    k = 0
    t_k = (E_start-E_rev)/scan
    print t_k
    while k < 2*l:
        if k < l:
            
            E_norm = 38.92*((E_start-E_o)-scan*k*(t_k)/(l))
            time.append((E_norm+E_o)/38.92)

            k_f = k_o*math.exp(-alpha*((E_norm-E_o)))
            k_b = k_o*math.exp((1-alpha)*((E_norm-E_o)))
            Z = (k_f*math.sqrt(t_k)/math.sqrt(Da)*CAold[0]-k_b*math.sqrt(t_k)/math.sqrt(Da)*CBold[0])
            if (Z*math.sqrt(DMA/(l))/2) > CAold[0]:
                current.append(CAold[0]*96485.3/t_k)
                CAnew.append(DMA*(CAold[1]-CAold[0]))
                if CBold[0]+DMB*(CBold[1]-CBold[0])+CAold[0]-k_rxn*CBold[0] > 0:
                    CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+CAold[0]-k_rxn*CBold[0])
                else:
                    CBnew.append(0)
            else:
                current.append(Z*math.sqrt(DMA/(l))*96485.3*math.sqrt(Da/t_k)/30)
                CAnew.append(CAold[0]+DMA*(CAold[1]-CAold[0])-Z*math.sqrt(DMA/(l)))
                if CBold[0]+DMB*(CBold[1]-CBold[0])+Z*math.sqrt(DMA/(l))-k_rxn*t_k/l*CBold[0] > 0:
                    CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+Z*math.sqrt(DMA/(l))-k_rxn*t_k/l*CBold[0])
                else:
                    CBnew.append(0)
            
            for i in range(1,int((l))):
                CAnew.append(CAold[i]+DMA*(CAold[i+1]-2*CAold[i]+CAold[i-1]))
                if CBold[i]+DMB*(CBold[i+1]-2*CBold[i]+CBold[i-1])-k_rxn*t_k/l*CBold[i] > 0:
                    CBnew.append(CBold[i]+DMB*(CBold[i+1]-2*CBold[i]+CBold[i-1])-k_rxn*t_k/l*CBold[i])
                else:
                    CBnew.append(0)

            CAnew.append(IniA)
            CBnew.append(IniB) 
            CAold = CAnew
            CAnew = []
            CBold = CBnew
            CBnew = []
            
            k = k+1
        else:
            E_norm = 38.92*((E_rev-E_o)+scan*(k-l)*(t_k)/(l))
            time.append((E_norm+E_o)/38.92)
            #Calculate current
            k_f = k_o*math.exp(-alpha*((E_norm-E_o)))
            k_b = k_o*math.exp((1-alpha)*((E_norm-E_o)))
            Z = (k_f*math.sqrt(t_k)/math.sqrt(Da)*CAold[0]-k_b*math.sqrt(t_k)/math.sqrt(Da)*CBold[0])
            if (Z*math.sqrt(DMA/(l))) > CAold[0]:
                current.append(CAold[0]*96485.3/t_k)
                CAnew.append(DMA*(CAold[1]-CAold[0]))
                CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+CAold[0]-k_rxn*t_k/l*CBold[0])
            else:
                current.append(Z*math.sqrt(DMA/(l))*96485.3*math.sqrt(Da/t_k)/30)
                CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+Z*math.sqrt(DMA/(l))-k_rxn*t_k/l*CBold[i])
                CAnew.append(CAold[0]+DMA*(CAold[1]-CAold[0])-Z*math.sqrt(DMA/(l)))
            
            #Diffusion beyond the first box
            for i in range(1,int((l))):
                CAnew.append(CAold[i]+DMA*(CAold[i+1]-2*CAold[i]+CAold[i-1]))
                CBnew.append(CBold[i]+DMB*(CBold[i+1]-2*CBold[i]+CBold[i-1])-k_rxn*t_k/l*CBold[i])

            CAnew.append(IniA)
            CBnew.append(IniB) 
            CAold = CAnew
            CAnew = []
            CBold = CBnew
            CBnew = []
            
            k = k+1
    master_time.append(time)
    master_current.append(current)

for i in range(len(master_current)):
    print max(master_current[i])
    print min(master_current[i])


#Current plot
fig = plt.figure()
ax = fig.add_subplot(111)
color_list = ['r','g','b','m','y']
for i in range(len(master_current)):
    ax.scatter(master_time[i],master_current[i],color=color_list[i],marker='o',s=1)

leg = ax.legend(('$k_{rxn} = 1 s^{-1}$','$k_{rxn} = 0.1 s^{-1}$','$k_{rxn} = 0.01 s^{-1}$','$k_{rxn} = 0.001 s^{-1}$','$k_{rxn} = 0.0001 s^{-1}$'), loc=2)
ax.set_xlim(ax.get_xlim()[::-1])
ax.set_xlabel('$E (V)$', fontsize="large")
ax.set_ylabel('$Current(arb.)$')
ax.legend()
plt.show()
