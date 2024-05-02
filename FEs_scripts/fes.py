#!/usr/bin/env python2
##
# Created: 13-Sep-2018 10:06:27 AM EDT
# Modified: 13-Sep-2018 01:54:02 PM EDT
# Created by: Matthew Varga
# Purpose: Plot 2D free energy surface from plumed sum_hills output
##

import matplotlib
import seaborn as sns
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import numpy as np
import sys

# Get data
print("Reading data...")
data = [line.split() for line in open(sys.argv[1], "r")]
data2 = [x for x in data if not x == []]  # strip out headers
file_name = str(sys.argv[1])
d1, d2, free, dd1, dd2 = [], [], [], [], []
for elem in data2[9:]:
    if elem[2] == 'inf' or elem[2] == '-inf':
        free.append(-10.0)
    else:
        free.append(float(elem[2]))

    d1.append(float(elem[0]))
    d2.append(float(elem[1]))
#    dd1.append(float(elem[3]))
#    dd2.append(float(elem[4]))

X = np.linspace(min(d1), max(d1), 1318)
Y = np.linspace(min(d2), max(d2), 1322)

print("Creating data grid. This may take a while...")
D1, D2 = np.meshgrid(X, Y)

#Normalize unbound state to zero
free = np.array(free)
# free = free+110.6

#Zero value
d1_arr = np.array(d1)
d2_arr = np.array(d2)
mask1 = (d1_arr >= 15.0 ) & (d1_arr < 16.0 )
mask2 = (d2_arr >= 15.0 ) & (d2_arr < 16.0 )
mask3 = mask1 * mask2

#correction = np.mean(free[mask3])
#free = free - correction

#Shift max value to a constant. To shifht the color map
#max_val = 50.0
#mask4 = (free>=max_val)
#free[mask4]=max_val

ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)
print(np.max(free))
print(np.min(free))
levels = np.arange(np.min(free), np.max(free), (np.max(free)-np.min(free))/20)
levels = levels.tolist()
levels = levels
levels = np.array(levels)
print(levels)
print(len(levels))
contour = plt.contour(D1, D2, ENER, colors='k', linewidths=0.3, levels=levels)
contourf = plt.contourf(D1, D2, ENER, cmap=cm.jet_r,levels=levels,alpha=0.98)
cbar = plt.colorbar()
plt.scatter([0.8525], [0.7766], marker='x', c='k')
plt.xlabel(r"$d_1$ (nm)",fontdict={'FontSize':20})
plt.ylabel(r"$d_2$ (nm)",fontdict={'FontSize':20})
plt.savefig(file_name+".png", dpi=300, bbox_inches='tight')
ax=plt.gca()
cbar.set_label(label="Free Energy (kJ/mol)",size='x-large',weight='bold')
ax.tick_params(axis='both',labelsize='xx-large')
cbar.ax.tick_params(axis='y',labelsize='x-large')
# sns.kdeplot(d1,d2,cmap='Spectral',shade=True)

ax.set_xlim(left=0.5,right=18.0)
ax.set_ylim(bottom=0.5,top=18.0)
plt.show()
"""
# 3D plots
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(D1, D2, ENER, cmap=cm.Spectral)
ax.set_zlim(np.min(ENER) - (np.max(ENER) - np.min(ENER))/10)
ax.contour(D1, D2, ENER, 10, linewidths=1, levels = levels, linestyles='solid', offset=ax.get_zlim()[0], cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='y', offset=np.min(D2), cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='x', offset=np.min(D1), cmap=cm.Spectral)
ax.set_zlabel("Free Energy (kJ/mol)")
ax.set_ylabel("Distance D2 (nm)")
ax.set_xlabel("Distance D1 (nm)")
ax.view_init(20, 45)
plt.draw()
plt.savefig("3dfes1.png", dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(D1, D2, ENER, cmap=cm.Spectral)
ax.set_zlim(np.min(ENER) - (np.max(ENER) - np.min(ENER))/10)
ax.contour(D1, D2, ENER, 10, linewidths=1, levels = levels, linestyles='solid', offset=ax.get_zlim()[0], cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='y', offset=np.max(D2), cmap=cm.Spectral)
ax.contour(D1, D2, ENER, zdir='x', offset=np.min(D1), cmap=cm.Spectral)
ax.set_zlabel("Free Energy (kJ/mol)")
ax.set_ylabel("Distance D2 (nm)")
ax.set_xlabel("Distance D1 (nm)")
ax.view_init(20, -45)
plt.draw()
plt.savefig("3dfes2.png", dpi=300, bbox_inches='tight')
"""
