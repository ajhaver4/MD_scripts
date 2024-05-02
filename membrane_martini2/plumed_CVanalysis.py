import numpy as np
import pandas
import matplotlib.pyplot as plt
import sys

data = pandas.read_csv(sys.argv[1],comment='#',delimiter='\s+',header=None,names=['time','d1','d2','minA','minB','membA','membB','l_angA','l_angB'])

fig,[ax,ax1,ax2] = plt.subplots(3,1)

ax.plot(data['time'],data['d1'],linewidth=0.6,color='royalblue',label='d1')
ax.plot(data['time'],data['d2'],linewidth=0.6,color='crimson',label='d2')
ax.set_xlabel("Time")
ax.set_ylabel("Distance (nm)")
ax.legend()

ax1.plot(data['time'],data['minA'],linewidth=0.6,color='orchid',label='Min dist A')
ax1.plot(data['time'],data['minB'],linewidth=0.6,color='peru',label='Min dist B')
ax1.set_xlabel("Time")
ax1.set_ylabel("Distance (nm)")
ax1.legend()

ax2.plot(data['time'],data['membA'],linewidth=0.6,color='orchid',label='Memb dist A')
ax2.plot(data['time'],data['membB'],linewidth=0.6,color='peru',label='Memb dist B')
ax2.set_xlabel("Time")
ax2.set_ylabel("Distance (nm)")
ax2.legend()

fig_ang,ax_ang = plt.subplots()
ax_ang.plot(data['time'],np.degrees(data['l_angA']),linewidth=0.6,color='forestgreen',label='l_ang A')
ax_ang.plot(data['time'],np.degrees(data['l_angB']),linewidth=0.6,color='gold',label='l_ang B')
ax_ang.set_xlabel("Time")
ax_ang.set_ylabel("Long Angle")
ax_ang.legend()

plt.show()
