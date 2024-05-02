import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

data = pd.read_csv(sys.argv[1],comment='#',delimiter="\s+",names=['time','d1','d2','min_A','min_B','lipid_A','lipid_B','langA','langB'])

n_frames = len(data['time'])
timesteps_array = data['time']

fig,ax = plt.subplots(3,1,sharex=True)
fig.suptitle("Distance Measurements from membrane")
ax[0].plot(timesteps_array,data['min_A'],label='Min Z dist A',color='forestgreen',linewidth=0.6)
ax[0].plot(timesteps_array,data['min_B'],label='Min Z dist B',color='cyan',linewidth=0.6)
ax[1].plot(timesteps_array,data['lipid_A'],label='Min dist A',color='forestgreen',linewidth=0.6)
ax[1].plot(timesteps_array,data['lipid_B'],label='Min dist B',color='cyan',linewidth=0.6)
ax[2].plot(timesteps_array,data['d1'],label='d1',color='blue',linewidth=0.6)
ax[2].plot(timesteps_array,data['d2'],label='d2',color='red',linewidth=0.6)

ax[2].set_xlabel("Time $\mu$s")
ax[0].set_ylabel("nm")
ax[1].set_ylabel("nm")
ax[2].set_ylabel("nm")
ax[0].legend()
ax[1].legend()
ax[2].legend()
# plt.ticklabel_format(axis='x',style='sci',scilimits=(0,3))

fig1,ax1 = plt.subplots()
ax1.plot(timesteps_array,180-np.degrees(data['langA']),label='l-angA',color='forestgreen',linewidth=0.6)
ax1.plot(timesteps_array,180-np.degrees(data['langB']),label='l-angB',color='cyan',linewidth=0.6)
ax1.set_xlabel("Time $\mu$s")
ax1.set_ylabel("deg")
ax1.legend()
plt.show()
