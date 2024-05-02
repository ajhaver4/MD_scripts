import sys
import numpy as np
import pandas # type: ignore
import matplotlib.pyplot as plt

data = pandas.read_csv(sys.argv[1],comment='#',sep='\t')
states_data = pandas.read_csv(sys.argv[2],delimiter='\s+')

#Creating states dictionary to hold boudnaries
states = {}
dG = {}
for i in range(len(states_data['states'])):

    coords = [float(states_data['xmin'][i]),float(states_data['xmax'][i]),float(states_data['ymin'][i]),float(states_data['ymax'][i])]
    states[states_data['states'][i]] = coords


n_values = len(data['d1'])


print("Selecting frames as per cut-off...")

cluster_time_dict = {}

fig,ax = plt.subplots()
for st,coords in states.items():
    x_min = coords[0]
    x_max = coords[1]
    y_min = coords[2]
    y_max = coords[3]
    # for i in range(0,n_values):
    #     if data['d1'][i] >= x_min and data['d1'][i] <= x_max and data['d2'][i] >= y_min and data['d2'][i] <= y_max:
    #         print(data['Timestep'][i],st)
    #         break

    mask1 = data['d1'] >= x_min
    mask2 = data['d1'] <= x_max
    mask3 = data['d2'] >= y_min
    mask4 = data['d2'] <= y_max

    filter = mask1 & mask2 & mask3 & mask4

    cluster_time_dict[st] = data['Timestep'][filter]

    ax.scatter(data['d1'][filter],data['d2'][filter],label=st)

ax.legend()

cry_timesteps=[]
for i in range(0,n_values):
    if data['RMSD-B'][i] < 1.5:
        cry_timesteps.append(data['Timestep'][i])

cluster_time_dict['CRY'] = cry_timesteps
ax.scatter(data['d1'][data['RMSD-B']<1.5],data['d2'][data['RMSD-B']<1.5],label='CRY')

plt.show()

print("Writing to file...")
with open("cluster_time_dict.dat",'w') as fl:
    for st,times in cluster_time_dict.items():
        fl.write(f"{st} {times}\n")



