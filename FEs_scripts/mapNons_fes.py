import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import numpy as np
import sys
import matplotlib.patches as patches
import pandas

# Get data
print("Reading data...")
data = [line.split() for line in open(sys.argv[1], "r")]
data2 = [x for x in data if not x == []]  # strip out headers
file_name = str(sys.argv[1])
states_data = pandas.read_csv(sys.argv[2],delimiter='\s+',comment='#')

d1, d2, free, dd1, dd2 = [], [], [], [], []
for elem in data2[9:]:
    d1.append(float(elem[0]))
    d2.append(float(elem[1]))
    free.append(float(elem[2]))
    # dd1.append(float(elem[3]))
    # dd2.append(float(elem[4]))

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

correction = np.mean(free[mask3])
free = free - correction

#Shift max value to a constant. To shifht the color map
max_val = 50.0
mask4 = (free>=max_val)
free[mask4]=max_val

ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)
print(np.max(free))
print(np.min(free))
levels = np.arange(np.min(free), np.max(free), (np.max(free)-np.min(free))/25)
levels = levels.tolist()
# levels = levels + [150.0]
levels = np.array(levels)
print(levels)
print(len(levels))
cmap = cm.get_cmap('jet_r')
cmap.set_over(color='black')
contour = plt.contour(D1, D2, ENER, colors='k', linewidths=0.3, levels=levels)
contourf = plt.contourf(D1, D2, ENER, cmap=cmap,levels=levels,alpha=0.98,vmax=max_val-15.0)
cbar = plt.colorbar()
plt.scatter([0.8525], [0.7766], marker='x', c='k',s=20)
plt.xlabel(r"$d_1$ (nm)",fontdict={'FontSize':20})
plt.ylabel(r"$d_2$ (nm)",fontdict={'FontSize':20})
plt.savefig(file_name+".png", dpi=300, bbox_inches='tight')
ax=plt.gca()
cbar.set_label(label="Free Energy (kJ/mol)",size=25,weight='bold')
ax.tick_params(axis='both',labelsize=20)
cbar.ax.tick_params(axis='y',labelsize=20)
#
#Format for each state box - [xmin,xmax,ymin,ymax]
#Endophilin
# states_cood = [[6.5,9.5,4.5,6.5],[13.0,15.0,4.0,6.0],[3.0,5.0,6.5,8.5],[1.5,2.5,4.0,6.5],[0.9,2.5,0.9,2.5]]
# states_cood = [[3.0,5.0,6.5,8.5],[0.9,2.5,0.9,2.5]]
#LSP1
# states_cood = [[2.5,4.5,0.8,1.3],[2.0,3.5,7.0,9.0],[6.5,8.0,3.0,4.0],[3.0,5.5,1.6,2.6],[0.9,2.5,0.9,2.5]]
# states_cood = [[2.0,3.5,7.0,9.0],[0.9,2.5,0.9,2.5]]
# label = ['NS1','NS2','NS3','NS4','Specific']

#Creating states dictionary to hold boudnaries
states = {}
dG = {}
for i in range(len(states_data['st'])):
    coords = [float(states_data['xmin'][i]),float(states_data['xmax'][i]),float(states_data['ymin'][i]),float(states_data['ymax'][i])]
    states[states_data['st'][i]] = coords
# label = ['NS3','Specific']
# color_code = {'NS1':'forestgreen','NS2':'crimson','NS3':'magenta','NS4':'orchid','Specific':'black'}

#color_st = {'CRY':'black','NS1':'magenta','NS5':'olive','TS1':'royalblue'}
#label_st = {'CRY':'Specific','NS1':'Non-Specific(1)','NS5':'Non-Specific(2)','TS1':'Two State'}

color_st = {'CRY':'black','NS1':'magenta','NS5':'olive','NS3':'indigo','TS2':'royalblue'}
label_st = {'CRY':'Specific','NS1':'Non-Specific(1)','NS5':'Non-Specific(2)','NS3':'Non-Specific(3)','TS2':'Two State'}

for st,cood in states.items():
    if ('unb' not in st) and ('TS1' not in st):
        wid = cood[1]-cood[0]
        ht = cood[3]-cood[2]
        rect = patches.Rectangle((cood[0],cood[2]),wid,ht,edgecolor=color_st[st],facecolor='None',linewidth=4.5,alpha=1.0)
        # ax.text(cood[0]+0.3,cood[2]+0.3,st,fontsize=12)
        ax.add_patch(rect)

# ax.hlines(12.5,0.8,12.5,colors='royalblue',linewidth=3.0)
# ax.vlines(12.5,0.8,12.5,colors='royalblue',linewidth=3.0)
ax.set_xlim(left=0.5,right=18.0)
ax.set_ylim(bottom=0.5,top=18.0)

val_range = np.arange(0,18,2)
ax.set_xticks(val_range)
ax.set_yticks(val_range)

ax.set_xticklabels(val_range,fontsize=20)
ax.set_yticklabels(val_range,fontsize=20)
# plt.savefig('States_boundary'+".png", dpi=300, bbox_inches='tight')
plt.show()
