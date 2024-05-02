import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")

n_frames = len(data['Timestep'])
timesteps_array = data['Timestep']/1e6
# timesteps_array = data['Frame']

fig,ax = plt.subplots(3,1,sharex=True)
fig.suptitle("Distance Measurements from membrane")
ax[0].plot(timesteps_array,data['min_dist'],label='Min dist',color='black',linewidth=0.6)
ax[0].plot(timesteps_array,data['min_dist_A'],label='Min dist chain A',color='forestgreen',linewidth=0.6)
ax[0].plot(timesteps_array,data['min_dist_B'],label='Min dist chain B',color='cyan',linewidth=0.6)
ax[1].plot(timesteps_array,data['d1'],label='d1',color='blue',linewidth=0.6)
ax[1].plot(timesteps_array,data['d2'],label='d2',color='red',linewidth=0.6)
#ax[2].plot(timesteps_array,data['membdist'],label='dist b/w MBS',color='gold',linewidth=0.6)
ax[2].plot(timesteps_array,data['memb_minA'],label='Memb_resid_A',color='darkorchid',linewidth=0.6)
ax[2].plot(timesteps_array,data['memb_minB'],label='Memb_resid_B',color='plum',linewidth=0.6)
ax[2].set_xlabel("Time $\mu$s")
ax[0].set_ylabel("nm")
ax[1].set_ylabel("nm")
ax[2].set_ylabel("nm")
ax[0].legend()
ax[1].legend()
ax[2].legend()
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,3))


fig2,ax2 = plt.subplots(2,1,sharex=True)
ax2[0].plot(timesteps_array,data['langle'],label='Long Angle',linewidth=0.6)
ax2[1].plot(timesteps_array,data['dim-angle'],label='Dim angle',linewidth=0.6)
ax2[1].set_xlabel("Time $\mu$s")
ax2[0].set_ylabel("degrees")
ax2[1].set_ylabel("degrees")
ax2[0].legend()
ax2[1].legend()
fig2.suptitle("Angle measurements")


fig3,ax3 = plt.subplots(2,1)
ax3[0].plot(timesteps_array,data['RMSD-A']/10.0,label="RMSD-A",linewidth=0.6,color='olive')
ax3[0].plot(timesteps_array,data['RMSD-B']/10.0,label='RMSD-B',linewidth=0.6,color='gold')
ax3[1].plot(timesteps_array,data['membdist']/10.0,label='Dist b/w memb patches',linewidth=0.6,color='orchid')
ax3[0].legend()
ax3[1].legend()
ax3[1].set_xlabel("Time $\mu$s")
ax3[0].set_ylabel("nm")
ax3[1].set_ylabel("nm")

# fig3,ax3 = plt.subplots(2,1)
# ax3[0].plot(timesteps_array,data['comA']/10.0,label="COM-A (z)",linewidth=0.6,color='olive')
# ax3[1].plot(timesteps_array,data['comB']/10.0,label='COM-B (z)',linewidth=0.6,color='purple')
# ax3[0].legend()
# ax3[1].legend()
# ax3[1].set_xlabel("Time $\mu$s")
# ax3[0].set_ylabel("nm")
# ax3[1].set_ylabel("nm")

#plt.show()

# comA_Z = []
# comB_Z = []
#
# for i in range(n_frames):
#     if timesteps_array[i] > 200000:
#         comA_Z.append(data['comA'][i]/10.0)
#         comB_Z.append(data['comB'][i]/10.0)

# #print(comA_Z)
# fig4,ax4 = plt.subplots()
# ax4.hist(comA_Z,bins=50,density=True,label="chainA")
# ax4.hist(comB_Z,bins=50,density=True,label="chainB")
# ax4.set_title("Distribution of COM (z) for chain A")
# ax4.legend()
# ax4.set_ylabel('P(z)')
# ax4.set_xlabel('Bins nm')

fig5,ax5 = plt.subplots()
ax5.plot(timesteps_array,data['zangle_A'],label="Chain A",linewidth=0.6,color='olive')
ax5.plot(timesteps_array,data['zangle_B'],label='Chain B',linewidth=0.6,color='purple')
ax5.legend()
ax5.set_xlabel("Time $\mu$s")
ax5.set_ylabel("degrees")
ax5.set_title("Angle between long axis and z-axis")
# fig,ax6 = plt.subplots(2,1,sharex=True)
# ax6[0].plot(timesteps_array,data['delta_z_A']/10.0,label='delta-chainA',color='black',linewidth=0.6)
# ax6[0].plot(timesteps_array,data['delta_z_B']/10.0,label='delta-chainB',color='blue',linewidth=0.6)
# ax6[1].plot(timesteps_array,data['lipid_count'],label='lipid_count',color='red',linewidth=0.6)
# ax6[1].set_xlabel("Time $\mu$s")
# ax6[0].set_ylabel("nm")
# ax6[1].set_ylabel("Number")
# ax6[0].legend()
# ax6[1].legend()

# fig7,ax7 = plt.subplots()
# ax7.plot(timesteps_array,data['wall_A']/10.0,label='Wall-chainA',color='black',linewidth=0.6)
# ax7.plot(timesteps_array,data['wall_B']/10.0,label='Wall-chainB',color='blue',linewidth=0.6)
# ax7.set_xlabel("Time $\mu$s")
# ax7.set_ylabel("nm")
# ax7.legend()
plt.show()
