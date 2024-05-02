import pandas
import matplotlib.pyplot as plt
import sys
import numpy as np

data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep='\s+',engine='python',names = ['time','d1','d2','bias','ct','rbias'])

time_array = np.arange(0,len(data['time']))*2
fig,ax = plt.subplots()

ax.plot(time_array[::10]/1e5,data['d1'][::10],color='red',label='d1',linewidth=0.6)
ax.plot(time_array[::10]/1e5,data['d2'][::10],color='blue',label='d2',linewidth=0.6)
ax.legend()
ax.set_xlabel(r'time $\mu$s')
ax.set_ylabel('nm')

# fig1,ax1 = plt.subplots()
# ax1.plot(time_array[:10000]/1e6,data['ct'][:10000],color='red',label='c(t)',linewidth=0.6)
# ax1.plot(time_array[:10000]/1e6,data['bias'][:10000],color='black',label='V(s,t)',linewidth=0.6)
# ax1.legend()
plt.show()
