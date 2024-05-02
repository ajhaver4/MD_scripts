import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import numpy as np
import sys
from scipy import constants
from matplotlib.widgets import Slider
# Get data
print("Reading data...")
data = [line.split() for line in open(sys.argv[1], "r")]
data2 = [x for x in data if not x == []]  # strip out headers
file_name = str(sys.argv[1])
d1, d2, free, dd1, dd2 = [], [], [], [], []
for elem in data2[9:]:
    d1.append(float(elem[0]))
    d2.append(float(elem[1]))
    free.append(float(elem[2]))
    dd1.append(float(elem[3]))
    dd2.append(float(elem[4]))
d1_bins = 879
d2_bins = 883


X = np.linspace(min(d1), max(d1), d1_bins)
Y = np.linspace(min(d2), max(d2), d2_bins)

#free_np = np.transpose(np.reshape(free,(d1_bins,d2_bins)))
free_np = np.reshape(free,(d2_bins,d1_bins))

print("Creating data grid. This may take a while...")
D1, D2 = np.meshgrid(X, Y)
#ENER = griddata((d1, d2), free, (D1, D2), method='linear', fill_value=0)

global major_d1
global major_d2

major_d1 = np.array(d1[:int(d1_bins)])
major_d2 = np.array(d2[:len(d2):int(d1_bins)])

def get_major_grid(rect,st):
    global major_d1
    global major_d2

    if st =='x':
        diff = major_d1 - rect
        bool = diff<=0
        bool_arr = diff[bool]
        #print(bool_arr)
        max = np.amax(bool_arr)
        idx = list(diff).index(max)
        return(idx)
    if st =='y':
        diff = major_d2 - rect
        bool = diff<=0
        bool_arr = diff[bool]
        max = np.amax(bool_arr)
        idx = list(diff).index(max)
        return(idx)

#Calculating the height of gaussina hills
h0 = 1     #Initial height - 1kJ/mol
biasF = 2    #(T + deltaT)/T = 10
#temp = 310
#deltaT = (biasF*temp) - temp
def calc_height(free_np,biasF,h0):
    temp = 310
    deltaT = (biasF*temp) - temp
    kb = constants.Boltzmann
    Na = constants.Avogadro
    h = h0*np.exp(free_np*1000/(kb*Na*deltaT))
    return(h)

h = calc_height(free_np,biasF,h0)

levels = np.arange(np.min(h), np.max(h), (np.max(h)-np.min(h))/36)
levels = levels.tolist()
#levels = levels
global contourf
global contour
global cb
fig,ax = plt.subplots()
contour = ax.contour(D1, D2, h, colors='k', linewidths=0.5, levels=levels)
contourf = ax.contourf(D1, D2, h, cmap=cm.Spectral, levels=levels)
cb = plt.colorbar(mappable=contourf,label="Height (kJ/mol)")
ax.set_xlabel("d1 (nm)")
ax.set_ylabel("d2 (nm)")
ax.set_title("Gaussian height")

#Make a slider
plt.subplots_adjust(bottom=0.25)
ax_slider = plt.axes([0.25,0.1,0.65,0.03])
bias_slider = Slider(ax=ax_slider,label='Bias factor',valmin=1,valmax=1000,valinit=biasF)

#Get heights for Two regions on FES
def get_height(h):
    ns_d1 = 2.2
    ns_d2 = 10.4
    sp_d1 = 0.872
    sp_d2 = 0.859

    ns_gridx = get_major_grid(ns_d1,'x')
    ns_gridy = get_major_grid(ns_d2,'y')
    sp_gridx = get_major_grid(sp_d1,'x')
    sp_gridy = get_major_grid(sp_d2,'y')

    ns_height = h[ns_gridy,ns_gridx]
    sp_height = h[sp_gridy,sp_gridx]

    return(ns_height,sp_height)

#Plot height difference and absolute height of crys
delta_H = []
cry_H = []
ns_H = []
bias_array = []
ns_h,sp_h = get_height(h)
delta_H.append(sp_h-ns_h)
bias_array.append(biasF)
cry_H.append(sp_h)
ns_H.append(ns_h)
def slider_update(val):
    global contourf
    global contour
    global cb
    for coll in contourf.collections:
        coll.remove()
    for coll in contour.collections:
        coll.remove()
    cb.remove()
    print("New bias factor: ",val)
    h = calc_height(free_np,val,h0)
    levels = np.arange(np.min(h), np.max(h), (np.max(h)-np.min(h))/36)
    levels = levels.tolist()
    contour = ax.contour(D1, D2, h, colors='k', linewidths=0.5, levels=levels)
    contourf = ax.contourf(D1, D2, h, cmap=cm.Spectral, levels=levels)
    cb = plt.colorbar(mappable=contourf,label="Height (kJ/mol)",ax=ax)

    plt.draw()
    print("Update heights")
    ns_h,sp_h = get_height(h)
    delta_H.append(sp_h-ns_h)
    bias_array.append(val)
    cry_H.append(sp_h)
    ns_H.append(ns_h)


bias_slider.on_changed(slider_update)


plt.show()

fig,ax = plt.subplots()
ax.plot(bias_array,delta_H,label='delta_H',linewidth=0.8)
ax.plot(bias_array,cry_H,label='Height at cry',linewidth=0.8)
ax.plot(bias_array,ns_H,label='Height at Non-specific',linewidth=0.8)
ax.axhline(1,linestyle='dashed',label=r'$\omega$',color='black')
ax.set_xlabel(r"Bias factor ($T + \Delta T$ / T)")
ax.set_ylabel("delta H (kJ/mol)")
ax.legend()
plt.show()
