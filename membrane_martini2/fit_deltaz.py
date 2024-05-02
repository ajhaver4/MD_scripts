import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")

n_frames = len(data['Timestep'])
timesteps_array = data['Timestep']/1e6
# timesteps_array = data['Frame']
threshold = data['Timestep'] > 0.1e5

fig,ax = plt.subplots()
ax.plot(timesteps_array,data['comA'],label='comA',linewidth=0.8,color='crimson')
ax.plot(timesteps_array,data['comB'],label='comB',linewidth=0.8,color='royalblue')
ax.set_ylabel("Distance (nm)")
ax.set_xlabel("Frame no. ")
ax.legend()

fig1,[ax1,ax2]=plt.subplots(2,1)
dzA=ax1.plot(timesteps_array,data['delta_z_A'],label='chain A',linewidth=0.8,color='crimson')
dzB=ax1.plot(timesteps_array,data['delta_z_B'],label='chain B',linewidth=0.8,color='royalblue')
fig1.suptitle("Delta Z values")
ax2.hist(data['delta_z_A'],bins=100,color='crimson',alpha=0.6,density=True)
ax2.hist(data['delta_z_B'],bins=100,color='royalblue',alpha=0.6,density=True)
ax1.set_ylabel("Distance (nm)")
ax1.set_xlabel("Frame no. ")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Distance (nm)")
ax1.legend()
ax2.legend()

fig3,[ax3,ax4,ax5]=plt.subplots(3,1)
ax3.hist(data['z_A_hi'],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax3.hist(data['z_B_hi'],bins=100,label='chain B',alpha=0.6,density=True,color='royalblue')
ax3.set_ylabel("Frequency")
ax3.set_xlabel("Distance (nm)")
ax3.set_title("High flexibility region")
ax3.legend()

ax4.hist(data['z_A_low'],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax4.hist(data['z_B_low'],bins=100,label='chain B',alpha=0.6,density=True,color='royalblue')
ax4.set_ylabel("Frequency")
ax4.set_xlabel("Distance (nm)")
ax4.set_title("Low flexibility region")
ax4.legend()

ax5.hist(data['z_A_med'],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax5.hist(data['z_B_med'],bins=100,label='chain B',alpha=0.6,density=True,color='royalblue')
ax5.set_ylabel("Frequency")
ax5.set_xlabel("Distance (nm)")
ax5.set_title("Mix flexibility region")
ax5.legend()

plt.show()


def harm2(x,k,r0):
    c = (2*3.14*2.5/k)**0.5
    v = 0.5*k*(x-r0)**2 + 2.5*np.log(c)
    return(v)


"""
Fitting membrane Distances
"""
figm,[axm1,axm2]=plt.subplots(2,1)
axm1.plot(timesteps_array[threshold],data['z_A_membres'][threshold],label='chain A',linewidth=0.8,color='crimson')
axm1.plot(timesteps_array[threshold],data['z_B_membres'][threshold],label='chain B',linewidth=0.8,color='steelblue')
fig1.suptitle("Membrane residues Delta Z values")

maskm = threshold & (data['z_A_membres'] < 1.3 ) & (data['z_A_membres'] > 0.7)  
memb_dist_A = data['z_A_membres'][maskm]

# maskb = (data['z_B_membres'] != 'nan') & threshold & (data['z_B_membres'] > 0.5 )
memb_dist_B = data['z_B_membres'][threshold]


pd_memb_A = axm2.hist(memb_dist_A,bins=100,color='crimson',alpha=0.6,density=True,label='ChainA')
pd_memb_B = axm2.hist(memb_dist_B,bins=100,color='steelblue',alpha=0.6,density=True,label='ChainB')


axm1.set_ylabel("Distance (nm)")
axm1.set_xlabel("Frame no. ")
axm2.set_ylabel("Frequency")
axm2.set_xlabel("Distance (nm)")
axm1.legend()
axm2.legend()

nzero_maskA = np.nonzero(pd_memb_A[0])
xdata_m = pd_memb_A[1][nzero_maskA]
ydata_m = pd_memb_A[0][nzero_maskA]


nzero_maskB = np.nonzero(pd_memb_B[0])
xdata_mB = pd_memb_B[1][nzero_maskB]
ydata_mB = pd_memb_B[0][nzero_maskB]


params_m = curve_fit(harm2,xdata_m,-2.5*np.log(ydata_m),[500,1.2])
params_mB = curve_fit(harm2,xdata_mB,-2.5*np.log(ydata_mB),[600,1.1])

figm3,axm3 = plt.subplots()
axm3.plot(xdata_m,-2.5*np.log(ydata_m),label='V(x)=-2.5 * log(P(x))')
axm3.plot(xdata_m,harm2(xdata_m,params_m[0][0],params_m[0][1]),label='V(x) = 0.5*k*x^2')
axm3.legend()
axm3.set_title("Chain A membrane Res")
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_m[0][0])
s2 = r'$mu$ = {:.2f}'.format(params_m[0][1])
axm3.text(np.mean(data['z_A_membres']),10,s1)
axm3.text(np.mean(data['z_A_membres']),11,s2)


figm4,axm4 = plt.subplots()
axm4.plot(xdata_mB,-2.5*np.log(ydata_mB),label='V(x)=-2.5 * log(P(x))')
axm4.plot(xdata_mB,harm2(xdata_mB,params_mB[0][0],params_mB[0][1]),label='V(x) = 0.5*k*x^2')
axm4.legend()
axm4.set_title("Chain B membrane Res")
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_mB[0][0])
s2 = r'$mu$ = {:.2f}'.format(params_mB[0][1])
axm4.text(np.mean(data['z_B_membres']),10,s1)
axm4.text(np.mean(data['z_B_membres']),11,s2)



k = np.linspace(10,200,10)
x=np.linspace(-2,2,50)


fig,ax = plt.subplots()
for i in range(len(k)):
    u=0.5*k[i]*(x**2)
    ax.plot(x,u,label=k[i],alpha=0.7,linewidth=0.8)

ax.legend()
ax.set_xlabel('r (nm)')
ax.set_ylabel('U(r)')
plt.show()
