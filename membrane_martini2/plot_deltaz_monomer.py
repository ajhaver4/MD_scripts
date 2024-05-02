import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

col = 'Frame\tTimestep\tcomA\tdelta_z_A\tz_A_hi\tz_A_low\tz_A_med\tz_A_membres'
col_names=col.split("\t")
print("Column Names: ",col_names)
# data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t",header=None,names=col_names)
data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")

n_frames = len(data['Timestep'])
timesteps_array = data['Timestep']/1e6

threshold = data['Timestep'] > 0.8e6
# timesteps_array = data['Frame']

fig,ax = plt.subplots()
ax.plot(timesteps_array,data['comA'],label='comA',linewidth=0.8,color='crimson')
ax.set_ylabel("Distance (nm)")
ax.set_xlabel("Frame no. ")
ax.legend()

fig1,[ax1,ax2]=plt.subplots(2,1)
ax1.plot(timesteps_array,data['delta_z_A'],label='chain A',linewidth=0.8,color='crimson')
fig1.suptitle("Delta Z values")
ax2.hist(data['delta_z_A'][threshold],bins=100,color='crimson',alpha=0.6,density=True,label='ChainA')
ax1.set_ylabel("Distance (nm)")
ax1.set_xlabel("Frame no. ")
ax2.set_ylabel("Frequency")
ax2.set_xlabel("Distance (nm)")
ax1.legend()
ax2.legend()

fig3,[ax3,ax4,ax5]=plt.subplots(3,1)
pd_hi_A=ax3.hist(data['z_A_hi'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax3.set_ylabel("Frequency")
ax3.set_xlabel("Distance (nm)")
ax3.set_title("High flexibility region")
ax3.legend()

pd_lo_A=ax4.hist(data['z_A_low'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax4.set_ylabel("Frequency")
ax4.set_xlabel("Distance (nm)")
ax4.set_title("Low flexibility region")
ax4.legend()
#
pd_mix_A=ax5.hist(data['z_A_med'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
ax5.set_ylabel("Frequency")
ax5.set_xlabel("Distance (nm)")
ax5.set_title("Mix flexibility region")
ax5.legend()

# plt.show()


#Curve fitting
def harm(x,k,r0):
    v = 0.5*k*(x-r0)**2
    return(v)

def harm2(x,k,r0):
    c = (2*3.14*2.5/k)**0.5
    v = 0.5*k*(x-r0)**2 + 2.5*np.log(c)
    return(v)

# print(len(pd_hi_A[1]))
# print(len(pd_hi_A[0]))
# sys.exit()


print("Average values: ")
mean_a = np.mean(data['z_A_low'])
print(mean_a)


mask1 = np.nonzero(pd_lo_A[0])
mask2 = np.nonzero(pd_hi_A[0])
mask3 = np.nonzero(pd_mix_A[0])
x_data = pd_lo_A[1][mask1]
y_data = pd_lo_A[0][mask1]
print(y_data)

params_01 = curve_fit(harm2,x_data,-2.5*np.log(y_data),[80,1.2])
params_02 = curve_fit(harm2,pd_hi_A[1][mask2],-2.5*np.log(pd_hi_A[0][mask2]),[12,1.0])
params_03 = curve_fit(harm2,pd_mix_A[1][mask3],-2.5*np.log(pd_mix_A[0][mask3]),[12,1.0])


fig_fit,[ax_fit2,ax_fit,ax_fit3] = plt.subplots(3,1)
ax_fit.plot(x_data,-2.5*np.log(y_data),label='V(x)=-2.5 * log(P(x))')
ax_fit.plot(x_data,harm2(x_data,params_01[0][0],mean_a),label='V(x) = 0.5*k*x^2')
ax_fit.legend()
ax_fit.set_title('Low flexibility Region')
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_01[0][0])
ax_fit.text(mean_a,10,s1)

ax_fit2.plot(pd_hi_A[1][mask2],-2.5*np.log(pd_hi_A[0][mask2]),label='V(x)=-2.5 * log(P(x))')
ax_fit2.plot(pd_hi_A[1][mask2],harm2(pd_hi_A[1][mask2],params_02[0][0],np.mean(data['z_A_hi'])),label='V(x) = 0.5*k*x^2')
ax_fit2.legend()
ax_fit2.set_title('High flexibility Region')
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_02[0][0])
ax_fit2.text(np.mean(data['z_A_hi']),10,s1)

ax_fit3.plot(pd_mix_A[1][mask3],-2.5*np.log(pd_mix_A[0][mask3]),label='V(x)=-2.5 * log(P(x))')
ax_fit3.plot(pd_mix_A[1][mask3],harm2(pd_mix_A[1][mask3],params_03[0][0],np.mean(data['z_A_med'])),label='V(x) = 0.5*k*x^2')
ax_fit3.legend()
ax_fit3.set_title('Mix flexibility Region')
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_03[0][0])
ax_fit3.text(np.mean(data['z_A_med']),10,s1)
fig_fit.suptitle("Chain A Fitting")


plt.show()


"""
Fitting membrane Distances
"""
figm,[axm1,axm2]=plt.subplots(2,1)
axm1.plot(timesteps_array,data['z_A_membres'],label='chain A',linewidth=0.8,color='crimson')
fig1.suptitle("Membrane residues Delta Z values")

maskm = (data['z_A_membres'] != 'nan') & threshold & (data['z_A_membres'] > 0.5 )
memb_dist_A = data['z_A_membres'][maskm]
pd_memb_A = axm2.hist(memb_dist_A,bins=100,color='crimson',alpha=0.6,density=True,label='ChainA')
axm1.set_ylabel("Distance (nm)")
axm1.set_xlabel("Frame no. ")
axm2.set_ylabel("Frequency")
axm2.set_xlabel("Distance (nm)")
axm1.legend()
axm2.legend()

nzero_maskA = np.nonzero(pd_memb_A[0])
xdata_m = pd_memb_A[1][nzero_maskA]
# ydata_m = pd_memb_A[0][nzero_maskA]*(pd_memb_A[1][1]-pd_memb_A[1][0])

# ydata_mB = pd_memb_B[0][nzero_maskB]*(pd_memb_B[1][1]-pd_memb_B[1][0])
ydata_m = pd_memb_A[0][nzero_maskA]

print(sum(pd_lo_A[0])*(pd_lo_A[1][1]-pd_lo_A[1][0]))
print(ydata_m)
# sys.exit()
params_m = curve_fit(harm2,xdata_m,-2.5*np.log(ydata_m),[300,1.0])

figm3,axm3 = plt.subplots()
axm3.plot(xdata_m,-2.5*np.log(ydata_m),label='V(x)=-2.5 * log(P(x))')
axm3.plot(xdata_m,harm2(xdata_m,params_m[0][0],params_m[0][1]),label='V(x) = 0.5*k*x^2')
axm3.legend()
axm3.set_title("Chain A membrane Res")
s1 = 'k = {:.2f} kJ mol-1 nm-2'.format(params_m[0][0])
s2 = r'$mu$ = {:.2f}'.format(params_m[0][1])
axm3.text(np.mean(data['z_A_membres']),10,s1)
axm3.text(np.mean(data['z_A_membres']),11,s2)



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
