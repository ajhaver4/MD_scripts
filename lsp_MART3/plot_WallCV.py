import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Plotting parameters
t_size=22
f_size=23
m_size=12
leg_size=24
fig_size=(6,5)

deltaZ_CV = pd.read_csv(sys.argv[1],comment='#',delimiter='\s+',names=['time','posA.z','posB.z','vol_res.bias','vol_res.force'])

fig,axA = plt.subplots(figsize=fig_size)
axA.plot(deltaZ_CV['time'],deltaZ_CV['posA.z'],label=r'$com_{A}$',color='blue')
axA.plot(deltaZ_CV['time'],deltaZ_CV['posB.z'],label=r'$com_{B}$',color='red')

axA.set_xlabel('Time (ps)',fontsize=f_size)
axA.set_ylabel('z (nm)',fontsize=f_size)
axA.legend(fontsize=leg_size)

fig,axB=plt.subplots(figsize=fig_size)
axB.plot(deltaZ_CV['time'],deltaZ_CV['vol_res.bias'],label=r'$V_{res}$',color='blue')

axB.set_xlabel('Time (ps)',fontsize=f_size)
axB.set_ylabel('Wall Pot (kJ/mol)',fontsize=f_size)
axB.legend(fontsize=leg_size)

plt.show()