from __future__ import division
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math


#Raw input setup
l = float(raw_input('Number of time steps.'))
DMA = float(raw_input('Dimensionless diffusion coefficient of species A.'))
DMB = float(raw_input('Dimensionless diffusion coefficient of species B.'))
Da = float(raw_input('Diffusion coefficient of species A.'))
IniA = float(raw_input('Molar concentration of species A.'))
IniB = float(raw_input('Molar concentration of species B.'))
k_o = float(raw_input('Electrode rate constant.'))
alpha = float(raw_input('Diffusion constant, alpha.'))
E_start = float(raw_input('Start potential.'))
E_rev = float(raw_input('End potential.'))
E_o = float(raw_input('Standard redox potential'))
scan = float(raw_input('Scan rate'))

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
while k < 2*l:
    if k < l:
        t_k = (E_start-E_rev)/scan
        E_norm = 38.92*((E_start-E_o)-scan*k*(t_k)/(l))
        time.append((E_norm+E_o)/38.92)

        k_f = k_o*math.exp(-alpha*((E_norm-E_o)))
        k_b = k_o*math.exp((1-alpha)*((E_norm-E_o)))
        Z = (k_f*math.sqrt(t_k)/math.sqrt(Da)*CAold[0]-k_b*math.sqrt(t_k)/math.sqrt(Da)*CBold[0])
        if (Z*math.sqrt(DMA/(l))) > CAold[0]:
            current.append(CAold[0])
            CAnew.append(DMA*(CAold[1]-CAold[0]))
            CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+CAold[0])
        else:
            current.append(Z*math.sqrt(DMA/(l)))
            CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+Z*math.sqrt(DMA/(l)))
            CAnew.append(CAold[0]+DMA*(CAold[1]-CAold[0])-Z*math.sqrt(DMA/(l)))
        
        for i in range(1,int((l))):
            CAnew.append(CAold[i]+DMA*(CAold[i+1]-2*CAold[i]+CAold[i-1]))
            CBnew.append(CBold[i]+DMB*(CBold[i+1]-2*CBold[i]+CBold[i-1]))

        CAnew.append(IniA)
        CBnew.append(IniB) 
        CAold = CAnew
        CAnew = []
        CBold = CBnew
        CBnew = []
        
        k = k+1
    else:
        t_k = (E_start-E_rev)/scan
        E_norm = 38.92*((E_rev-E_o)+scan*(k-l)*(t_k)/(l))
        time.append((E_norm+E_o)/38.92)
        #Calculate current
        k_f = k_o*math.exp(-alpha*((E_norm-E_o)))
        k_b = k_o*math.exp((1-alpha)*((E_norm-E_o)))
        Z = (k_f*math.sqrt(t_k)/math.sqrt(Da)*CAold[0]-k_b*math.sqrt(t_k)/math.sqrt(Da)*CBold[0])
        if (Z*math.sqrt(DMA/(l))) > CAold[0]:
            current.append(CAold[0])
            CAnew.append(DMA*(CAold[1]-CAold[0]))
            CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+CAold[0])
        else:
            current.append(Z*math.sqrt(DMA/(l)))
            CBnew.append(CBold[0]+DMB*(CBold[1]-CBold[0])+Z*math.sqrt(DMA/(l)))
            CAnew.append(CAold[0]+DMA*(CAold[1]-CAold[0])-Z*math.sqrt(DMA/(l)))
        
        #Diffusion beyond the first box
        for i in range(1,int((l))):
            CAnew.append(CAold[i]+DMA*(CAold[i+1]-2*CAold[i]+CAold[i-1]))
            CBnew.append(CBold[i]+DMB*(CBold[i+1]-2*CBold[i]+CBold[i-1]))

        CAnew.append(IniA)
        CBnew.append(IniB) 
        CAold = CAnew
        CAnew = []
        CBold = CBnew
        CBnew = []
        
        k = k+1


print max(current)
print min(current)
#Current plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(time,current,c='r',marker='^')
ax.set_xlim(ax.get_xlim()[::-1])
ax.set_xlabel('$E (V)$', fontsize="large")
ax.set_ylabel('$Current(arb.)$')
ax.legend()
plt.show()
