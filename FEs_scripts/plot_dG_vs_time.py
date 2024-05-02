import pandas
import sys
import numpy as np
import matplotlib.pyplot as plt

data = pandas.read_csv(sys.argv[1],comment='#',delimiter='\t')
states_data = pandas.read_csv(sys.argv[2],delimiter='\s+')
states={}
for i in range(len(states_data['st'])):

    states[states_data['st'][i]] = 0.0

#Membrane
V = 37.88*37.88*18.03   #nm3
#V = 37.88*37.88*37.83   #nm3
c0 = 1e3 * (1/1660)    #nm3

std_c = -8.314*310*np.log(V/c0)/1000    #kJ/mol

print("Standard State Correction :", std_c)
# sys.exit()
fig,ax=plt.subplots()

color_st = {'CRY':'black','NS1':'magenta','NS2':'olive','NS3':'steelblue','TS1':'royalblue'}
label_st = {'CRY':'Specific','NS1':'Non-Specific(1)','NS2':'Non-Specific(2)','NS3':'Non-Specific(3)','TS1':r'$\Delta G_{2D} $'}

for st,v in states.items():
    if ('unb' not in st) and ('TS2' not in st):
    #if ('TS2' in st):
        # print(st,data[st])
        new_data = data[st] > -500
        ax.plot(data['Timestep']/1e6,data[st],linewidth=2.5,alpha=0.8,label=label_st[st],color=color_st[st])
        # ax.plot(data['Timestep'][new_data]/1e6,data[st][new_data],linewidth=4.0,alpha=0.8,label=label_st[st],color=color_st[st])
        print("State: ", st)
        print("Free Energy : ",data[st].to_numpy()[-1])

        kd = 1e9*np.exp((data[st].to_numpy()[-1])*1000/(8.314*310))
        #kd = 1e9*np.exp((data[st].to_numpy()[-1])*1000/(8.314*310))
        # kd=0
        print("State: ", st)
        print("Binding Affinity: %.2f nM \t\t dG: %.2f" %(kd,data[st].to_numpy()[-1]))

ax.legend(loc='best',fontsize=20)
ax.set_ylabel(r'$\Delta G$ kJ/mol',fontsize=30)
ax.set_xlabel(r"Time in $\mu s$",fontsize=30)
ax.tick_params(labelsize=20)
plt.show()
