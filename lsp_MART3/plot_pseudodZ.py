import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

z_data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")
final_data = pd.read_csv(sys.argv[2],comment='#',sep='\t')


#Plotting parameters
t_size=22
f_size=23
m_size=12
leg_size=24
fig_size=(6,5)

#BOUND STATE FILTERING
filter_time = []
for i in range(len(final_data['d1'])):
    if final_data['d1'][i]<2.0 or final_data['d1'][i]<2.0:
        filter_time.append(final_data['Timestep'][i])

dataA = z_data[z_data['Timestep'].isin(filter_time)]
dataB = z_data[z_data['Timestep'].isin(filter_time)]

thresholdA = dataA['Timestep'] > 0.5e6
thresholdB = dataB['Timestep'] > 0.5e6

deltaZ_CV = pd.read_csv(sys.argv[3],comment='#',delimiter='\s+',names=['time','p1.z','p2.z','r1.bias','r2.bias'])

fig,axA = plt.subplots(figsize=fig_size)



p4A=axA.hist(dataA['z_A_membres'][thresholdA]-np.mean(dataA['z_A_membres'][thresholdA]),bins=100,label=r'$P_{mem} (\psi)$',alpha=0.6,density=True,color='skyblue')
p4A1 = axA.hist(deltaZ_CV['p1.z'] - np.mean(deltaZ_CV['p1.z']),bins=100,label=r'$P_{res} (\psi)$',alpha=0.4,density=True,color='k')
fig.tight_layout()
axA.tick_params(labelsize=t_size)
axA.set_xlabel(r'$dz$ (nm)',fontsize=f_size)
axA.set_ylabel(r"$P (dz)$",fontsize=f_size)
axA.legend(fontsize=leg_size)

fig,axB = plt.subplots(figsize=fig_size)
p4B=axB.hist(dataA['z_B_membres'][thresholdB] - np.mean(dataA['z_B_membres'][thresholdB]),bins=100,label=r'$P_{mem} (\psi)$',alpha=0.6,density=True,color='skyblue')
p4B1 = axB.hist(deltaZ_CV['p2.z'] - np.mean(deltaZ_CV['p2.z']),bins=100,label=r'$P_{res} (\psi)$',alpha=0.4,density=True,color='k')
fig.tight_layout()
axB.tick_params(labelsize=t_size)
axB.set_xlabel(r'$dz$ (nm)',fontsize=f_size)
axB.set_ylabel(r"$P (dz)$",fontsize=f_size)
axB.legend(fontsize=leg_size)

plt.show()

