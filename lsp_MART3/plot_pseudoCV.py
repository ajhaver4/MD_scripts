import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

pitch_data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")
final_data = pd.read_csv(sys.argv[2],comment='#',sep='\t')


#Plotting parameters
t_size=22
f_size=23
m_size=12
leg_size=24
fig_size=(6,5)

#For BOUND STATE FILTERING
filter_time = []

for i in range(len(final_data['d1'])):
    if final_data['d1'][i]<2.0 or final_data['d1'][i]<2.0:
        filter_time.append(final_data['Timestep'][i])
# data=pitch_data[pitch_data['Timestep'].isin(filter_time)]
dataA=pitch_data[pitch_data['Timestep'].isin(filter_time)]
dataB=pitch_data[pitch_data['Timestep'].isin(filter_time)]

thresholdA = dataA['Timestep'] > 0.5e6
thresholdB = dataB['Timestep'] > 0.5e6

pseudo_CV = pd.read_csv(sys.argv[3],comment='#',delimiter='\s+',names=['time','pitchA','pitchB','rollA','rollB','res_pitchA.bias', 'res_pitchB.bias', 'res_rollA.bias','res_rollB.bias'])


fig,axA = plt.subplots(figsize=fig_size)

p4A=axA.hist(dataA['PitchA3'][thresholdA],bins=100,label=r'$P_{mem} (\psi)$',alpha=0.6,density=True,color='limegreen')
p4A1 = axA.hist(np.degrees(pseudo_CV['pitchA']),bins=100,label=r'$P_{res} (\psi)$',alpha=0.4,density=True,color='k')
fig.tight_layout()
axA.tick_params(labelsize=t_size)
axA.set_xlabel(r'$\psi$ (deg)',fontsize=f_size)
axA.set_ylabel(r"$P_{mem} (\psi)$",fontsize=f_size)
axA.legend(fontsize=leg_size)

fig,axB = plt.subplots(figsize=fig_size)
p4B=axB.hist(dataB['PitchB3'][thresholdB],bins=100,label=r'$P_{mem} (\psi)$',alpha=0.6,density=True,color='limegreen')
p4B1 = axB.hist(np.degrees(pseudo_CV['pitchB']),bins=100,label=r'$P_{res} (\psi)$',alpha=0.4,density=True,color='k')
fig.tight_layout()
axB.tick_params(labelsize=t_size)
axB.set_xlabel(r'$\psi$ (deg)',fontsize=f_size)
axB.set_ylabel(r"$P_{mem} (\psi)$",fontsize=f_size)
axB.legend(fontsize=leg_size)

fig,axA = plt.subplots(figsize=fig_size)
p4A=axA.hist(dataA['RollA4'][thresholdA],bins=100,label=r'$P_{mem} (\phi)$',alpha=0.6,density=True,color='darkorange')
p4A2 = axA.hist(np.degrees(pseudo_CV['rollA']),bins=100,label=r'$P_{res} (\phi)$',alpha=0.4,density=True,color='k')
fig.tight_layout()
axA.tick_params(labelsize=t_size)
axA.set_xlabel(r'$\phi$ (deg)',fontsize=f_size)
axA.set_ylabel(r"$P_{mem} (\phi)$",fontsize=f_size)
axA.legend(fontsize=leg_size)

fig,axB = plt.subplots(figsize=fig_size)
p4B=axB.hist(dataB['RollB4'][thresholdB],bins=100,label=r'$P_{mem} (\phi)$',alpha=0.6,density=True,color='darkorange')
p4B2 = axB.hist(np.degrees(pseudo_CV['rollB']),bins=100,label=r'$P_{res} (\phi)$',alpha=0.4,density=True,color='k')
fig.tight_layout()
axB.tick_params(labelsize=t_size)
axB.set_xlabel(r'$\phi$ (deg)',fontsize=f_size)
axB.set_ylabel(r"$P(\phi)$",fontsize=f_size)
axB.legend(fontsize=leg_size)

plt.show()


