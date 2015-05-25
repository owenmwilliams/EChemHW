from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import math



scan_list = [0.1,1,10]


#Array setup          
dist=[]
conc_final = []
master_pot = []
master_current = []
time = []
time_list = []
#Iterative Calculations for Current/Concentration

for i in range(len(scan_list)):
    l = 100
    Da = float(0.00001)
    phi = float(0.0001)
    DMA = float(0.45)
    E_start = float(0.25)
    E_rev = float(-0.25)
    area = float(1)
    alpha = float(0.5)
    vol = float(phi*area)
    scan = float(scan_list[i])
    t = float((E_start-E_rev)/scan)
    delta_t = float(t/l)
    delta_x = math.sqrt(Da*delta_t/DMA)
    
    time_list.append(t)
    time_list
    if int(delta_x/phi) >= 1:
        print 'yes'
        
        conc = float(0.000001)
        conc2 = 0
        potential = []
        current = []
        current2 = []
        potential2 = []
        k = 0
        while k < 100:
        
            dist.append(delta_x*k)
            E = E_start - delta_t*k*scan
            potential.append(E)
            potential2.append(E)
            current_o = 96485.3*38.92*scan*vol*conc*(math.exp(38.92*E)/(1+math.exp(38.92*E))**2)
            if current_o*delta_t/(96485.3*vol) > conc:
                current.append(conc*vol*96485.3/delta_t)
                current2.append(-conc*vol*96485.3/delta_t)


                conc2 = conc2 + conc
                conc = 0
            else:
                current.append(current_o)
                current2.append(-current_o)

                conc = conc - current_o*delta_t/(96485.3*vol)
                conc2 = conc2 + current_o*delta_t/(96485.3*vol)

            k = k+1
        master_pot.append(potential)
        master_pot.append(potential2)
    
        master_current.append(current)
        master_current.append(current2)

    else:
        l = int(phi/delta_x)
        if l == 1:
            print 'a'
            l = 2
            delta_t = float(t/2)
            delta_x = math.sqrt(Da*delta_t/DMA)
            print delta_x
            CAold = []
            CAnew = []
            CBold = []
            CBnew = []
            conc = float(0.000001)
            for i in range(l):
                CAold.append(conc)
                CBold.append(0)
            potential = []
            current = []
            k = 0
            while k < l:
                E = E_start - delta_t*k*scan
                potential.append(E)
                current_o = 96485.3*38.92*scan*area*delta_x*CAold[0]*(math.exp(38.92*E)/(1+math.exp(38.92*E))**2)
                if current_o*delta_t/(96485.3*area*delta_x) < DMA*(CAold[1]-CAold[0]):
                    current.append(current_o)
                    CAnew.append(DMA*(CAold[1]-CAold[0]) - current_o*delta_t/(96485.3*area*delta_x))
                    CAnew.append(CAold[1]+DMA*(CAold[0]-CAold[1]))
                    CBnew.append(DMA*(CBold[1]-CBold[0]) + current_o*delta_t/(96485.3*area*delta_x))
                    CBnew.append(CBold[1]+DMA*(CBold[0]-CBold[1]))
                else:
                    current.append(CAold[0]*area*delta_x*96485.3/delta_t)
                    CAnew.append(DMA*(CAold[1]-CAold[0]))
                    CAnew.append(CAold[1] - DMA*(CAold[1]-CAold[0]))
                    CBnew.append(DMA*(CBold[1]-CBold[0]) + CAold[0]*area*delta_x*96485.3/delta_t)
                    CBnew.append(CBold[1]+DMA*(CBold[0]-CBold[1]))
                CAold = CAnew
                CAnew = []
                CBold = CBnew
                CBnew = []
                k = k+1
            master_pot.append(potential)
            master_current.append(current)

        else:
            print 'b'
            delta_t = float(t/l)
            delta_x = math.sqrt(Da*delta_t/DMA)
            print delta_x
            CAold = []
            CAnew = []
            CBold = []
            CBnew = []
            conc = float(0.000001)
            for i in range(l):
                CAold.append(conc)
                CBold.append(0)
            potential = []
            current = []
            k = 0
            while k < l:
                E = E_start - delta_t*k*scan
                potential.append(E)
                current_o = 96485.3*38.92*scan*area*delta_x*CAold[0]*(math.exp(38.92*E)/(1+math.exp(38.92*E))**2)
                if current_o*delta_t/(96485.3*area*delta_x) < DMA*(CAold[1]-CAold[0]):
                    current.append(current_o)
                    CAnew.append(DMA*(CAold[1]-CAold[0]) - current_o*delta_t/(96485.3*area*delta_x))
                    CBnew.append(DMA*(CBold[1]-CBold[0]) + current_o*delta_t/(96485.3*area*delta_x))
                    for i in range(l-2):
                        CAnew.append(CAold[i+1] + DMA*(CAold[i]-2*CAold[i+1]+CAold[i+2]))
                        CBnew.append(CBold[i+1] + DMA*(CBold[i]-2*CBold[i+1]+CBold[i+2]))
                    CBnew.append(CBold[1]+DMA*(CBold[0]-CBold[1]))
                    CAnew.append(CAold[1]+DMA*(CAold[0]-CAold[1]))
                else:
                    current.append(CAold[0]*area*delta_x*96485.3/delta_t)
                    CAnew.append(DMA*(CAold[1]-CAold[0]))
                    CBnew.append(DMA*(CBold[1]-CBold[0]) + CAold[0]*area*delta_x*96485.3/delta_t)
                    for i in range(l-2):
                        CAnew.append(CAold[i+1] + DMA*(CAold[i]-2*CAold[i+1]+CAold[i+2]))
                        CBnew.append(CBold[i+1] + DMA*(CBold[i]-2*CBold[i+1]+CBold[i+2]))
                    CAnew.append(CAold[1] - DMA*(CAold[1]-CAold[0]))
                    CBnew.append(CBold[1]+DMA*(CBold[0]-CBold[1]))
                CAold = CAnew
                CAnew = []
                CBold = CBnew
                CBnew = []
                k = k+1
            master_pot.append(potential)
            master_current.append(current)
            
  
#Current plot
fig = plt.figure()
ax = fig.add_subplot(111)
color_list = ['r','g','b','m']
for i in range(len(master_current)):
    ax.plot(master_pot[i],master_current[i],color = color_list[int(i/2)])


ax.set_xlim(ax.get_xlim()[::-1])
ax.set_xlabel('$E (V)$', fontsize="large")
ax.set_ylabel('$Current$')
ax.legend()
plt.grid(True)
plt.show()

