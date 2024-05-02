import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

pitch_data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")
final_data = pd.read_csv(sys.argv[2],comment='#',sep='\t')

t_size=22
f_size=23
m_size=12
leg_size=24
fig_size=(6,5)

#For BOUND STATE FILTERING
filter_time = []

for i in range(len(final_data['d1'])):
    if final_data['d1'][i]<2.0 or final_data['d1'][i]<2.0:
        filter_time.append(final_data['Timestep'][i])
# data=pitch_data[pitch_data['Timestep'].isin(filter_time)]
dataA=pitch_data[pitch_data['Timestep'].isin(filter_time)]
dataB=pitch_data[pitch_data['Timestep'].isin(filter_time)]


#FOR UNBOUND STATE FILTERING
# filter_time_A=[]
# filter_time_B=[]
# for i in range(len(final_data['d1'])):
#     if final_data['d1'][i]>12.0 and final_data['d1'][i]>12.0:
#         if final_data['min_dist_A'][i] <= 1.0:
#             filter_time_A.append(final_data['Timestep'][i])
#         if final_data['min_dist_B'][i] <=1.0:
#             filter_time_B.append(final_data['Timestep'][i])

# dataA = pitch_data[pitch_data['Timestep'].isin(filter_time_A)]
# dataB = pitch_data[pitch_data['Timestep'].isin(filter_time_B)]


# n_frames = len(data['Timestep'])
timesteps_arrayA = dataA['Timestep']/1e6
timesteps_arrayB = dataB['Timestep']/1e6
thresholdA = dataA['Timestep'] > 0.5e6
thresholdB = dataB['Timestep'] > 0.5e6

figp,[axp1,axp2,axp3,axp4,axp5] = plt.subplots(5,1)
p1A = axp1.hist(dataA['Z_angle_A'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p1B = axp1.hist(dataB['Z_angle_B'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp1.set_title("Z-angle")
axp1.legend()

p2A=axp2.hist(dataA['PitchA1'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p2B=axp2.hist(dataB['PitchB1'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp2.set_title("Pitch 1")
axp2.legend()
p3A=axp3.hist(dataA['PitchA2'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p3B=axp3.hist(dataB['PitchB2'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp3.set_title("Pitch 2")
axp3.legend()

p4A=axp4.hist(dataA['PitchA3'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p4B=axp4.hist(dataB['PitchB3'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp4.set_title("Pitch 3")
axp4.legend()

p5A=axp5.hist(dataA['PitchA4'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p5B=axp5.hist(dataB['PitchB4'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp5.set_title("Pitch 4")
axp5.legend()

figr,[axr1,axr2,axr3,axr4,axr5] = plt.subplots(5,1)
r1A = axr1.hist(dataA['RollA1'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r1B = axr1.hist(dataB['RollB1'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr1.set_title("Roll 1")
axr1.legend()

r2A=axr2.hist(dataA['RollA2'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r2B=axr2.hist(dataB['RollB2'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr2.set_title("Roll 2")
axr2.legend()

r3A=axr3.hist(dataA['RollA3'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r3B=axr3.hist(dataB['RollB3'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr3.set_title("Roll 3")
axr3.legend()

r4A=axr4.hist(dataA['RollA4'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r4B=axr4.hist(dataB['RollB4'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr4.set_title("Roll 4")
axr4.legend()

r5A=axr5.hist(dataA['RollA5'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r5B=axr5.hist(dataB['RollB5'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr5.set_title("Roll 5")
axr5.legend()


figA,[axA,axA1] = plt.subplots(2,1,sharex=True)
axA.plot(timesteps_arrayA,dataA['Z_angle_A'],label='Z-angle',linewidth=0.7)
axA.plot(timesteps_arrayA,dataA['PitchA1'],label='PitchA1',linewidth=0.7)
axA.plot(timesteps_arrayA,dataA['PitchA2'],label='PitchA2',linewidth=0.7)
axA.plot(timesteps_arrayA,dataA['PitchA3'],label='PitchA3',linewidth=0.7)
axA.plot(timesteps_arrayA,dataA['PitchA4'],label='PitchA4',linewidth=0.7)
axA1.plot(timesteps_arrayB,dataB['Z_angle_B'],label='Z-angle',linewidth=0.7)
axA1.plot(timesteps_arrayB,dataB['PitchB1'],label='PitchB1',linewidth=0.7)
axA1.plot(timesteps_arrayB,dataB['PitchB2'],label='PitchB2',linewidth=0.7)
axA1.plot(timesteps_arrayB,dataB['PitchB3'],label='PitchB3',linewidth=0.7)
axA1.plot(timesteps_arrayB,dataB['PitchB4'],label='PitchB4',linewidth=0.7)
axA.legend()
axA1.legend()
axA.set_xlabel("Time")
axA.set_ylabel("Angle ")
axA1.set_ylabel("Angle ")
axA.set_title("Chain A")
axA1.set_title("Chain B")

figB,[axB,axB1] = plt.subplots(2,1,sharex=True)
axB.plot(timesteps_arrayA,dataA['RollA1'],label="Roll1",linewidth=0.7)
axB.plot(timesteps_arrayA,dataA['RollA2'],label="Roll2",linewidth=0.7)
axB.plot(timesteps_arrayA,dataA['RollA3'],label="Roll3",linewidth=0.7)
axB.plot(timesteps_arrayA,dataA['RollA4'],label="Roll4",linewidth=0.7)
axB.plot(timesteps_arrayA,dataA['RollA5'],label="Roll5",linewidth=0.7)
axB1.plot(timesteps_arrayB,dataB['RollB1'],label="Roll1",linewidth=0.7)
axB1.plot(timesteps_arrayB,dataB['RollB2'],label="Roll2",linewidth=0.7)
axB1.plot(timesteps_arrayB,dataB['RollB3'],label="Roll3",linewidth=0.7)
axB1.plot(timesteps_arrayB,dataB['RollB4'],label="Roll4",linewidth=0.7)
axB1.plot(timesteps_arrayB,dataB['RollB5'],label="Roll5",linewidth=0.7)
axB.legend()
axB1.legend()
axB.set_xlabel("Time")
axB.set_ylabel("Angle ")
axB1.set_ylabel("Angle ")
axB.set_title("Chain A")
axB1.set_title("Chain B")
plt.show()


#FItting
def harm2(theta,k):
    v = k*(1-np.cos(1*(np.radians(theta))))
    return(v)

def harm(theta,k,theta0):
    c = (2*3.14*2.5/k)**0.5
    v = 0.5*k*(np.radians(theta-theta0))**2 #+ 2.5*np.log(c)
    return(v)

"""Chain A Pitch"""

hist_pitchA = p4A

mask1 = np.nonzero(hist_pitchA[0])
avg_pitch = np.mean(hist_pitchA[1][mask1])
# mask1 = (hist_pitchA[1] >= 65) & (hist_pitchA[1] <=100)
# mask1 = hist_pitchA[0]>0.002


# x_data = p4A[1][mask1]-avg_pitch
x_data = hist_pitchA[1][mask1]
y_data = -2.5*np.log(hist_pitchA[0][mask1])+2.5*np.log(np.max(hist_pitchA[0][mask1]))
# y_data = -2.5*np.log(p4A[0][mask1])
params_02 = curve_fit(harm,x_data,y_data,[400,75],maxfev=8000,bounds=(1.0,2000))
#
#
#
fig_fit,ax_fit = plt.subplots(figsize=fig_size)
ax_fit.plot(x_data,y_data,label=r'$U_{mem} (\psi)$',alpha=0.7,color='limegreen',linestyle='dashed',marker='o',markersize=m_size)
ax_fit.plot(x_data,harm(x_data,*params_02[0]),label=r'$U_{fit} (\psi)$',alpha=0.7,color='k')
ax_fit.legend(fontsize=leg_size)
# ax_fit.set_title(r'Pitch ($\psi$) - Chain A')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02[0][0])
# s2 = 'n = {:.2f}'.format(params_02[0][1])
s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
# ax_fit.text(np.mean(x_data),10,s1)
# ax_fit.text(np.mean(x_data),12,s2)
# ax_fit.text(0,14,s2)
fig_fit.tight_layout()
ax_fit.tick_params(labelsize=t_size)
ax_fit.set_xlabel(r'$\psi$ (deg)',fontsize=f_size)
ax_fit.set_ylabel("Potential Energy (kJ/mol)",fontsize=f_size)

"""Chain B Pitch"""

hist_pitchB = p4B

mask2 = np.nonzero(hist_pitchB[0])

# mask2 = (hist_pitchB[1] >= 52) & (hist_pitchB[1] <=106)
avg_pitchB = np.mean(p4B[1][mask2])
# x_data1 = p4B[1][mask2]-avg_pitchB
x_data1 = hist_pitchB[1][mask2]
y_data1 = -2.5*np.log(hist_pitchB[0][mask2]) + 2.5*np.log(np.max(hist_pitchB[0][mask2]))
# y_data1 = -2.5*np.log(p4B[0][mask2])
# print("Pitch Shift U(mean): ", 2.5*np.log(np.max(hist_pitchB[0][mask2])))
params_02B = curve_fit(harm,x_data1,y_data1,[450,10],maxfev=8000,bounds=(1.0,2000))

print("Pitch Data")
print("Avg for chainA : ",np.radians(avg_pitch))
print("Avg for chainB: ", np.radians(avg_pitchB))

fig_fit,ax_fit2 = plt.subplots(figsize=fig_size)
ax_fit2.plot(x_data1,y_data1,label=r'$U_{mem} (\psi)$',alpha=0.7,color='limegreen',linestyle='dashed',marker='o',markersize=m_size)
# ax_fit2.plot(x_data1,harm(x_data1,params_02[0][0]),label='V(x) = k*(1-cos(theta)) for A')
ax_fit2.plot(x_data1,harm(x_data1,*params_02B[0]),label=r'$U_{fit} (\psi)$',alpha=0.7,color='k')
ax_fit2.legend(fontsize=leg_size)
# ax_fit2.set_title(r'Pitch ($\psi$) - Chain B')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02B[0][0])
# s2 = 'n = {:f}'.format(params_02B[0][1])
s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitchB))
# ax_fit2.text(np.mean(x_data1),10,s1)
# ax_fit2.text(np.mean(x_data1),12,s2)
# ax_fit2.text(0,14,s2)
fig_fit.tight_layout()
ax_fit2.tick_params(labelsize=t_size)
ax_fit2.set_xlabel(r'$\psi$ (deg)',fontsize=f_size)
ax_fit2.set_ylabel("Potential Energy (kJ/mol)",fontsize=f_size)

"""Chain A Roll"""
mask1 = np.nonzero(r4A[0])
# mask1 = (r3A[1] >=30) & (r3A[1] <=91)
avg_roll = np.mean(r4A[1][mask1])
# x_data = r4A[1][mask1]-avg_pitch
x_data = r4A[1][mask1]
y_data = -2.5*np.log(r4A[0][mask1])+2.5*np.log(np.max(r4A[0][mask1]))
# y_data = -2.5*np.log(r4A[0][mask1])
params_02A = curve_fit(harm,x_data,y_data,[230,18])

"""Chain B Roll"""
mask2 = r4B[0] > 0
print("Length of Hist otput: ",len(r4B[0]),len(r4B[1]))
mask3 = (r4B[1][1:] >2) & (r4B[1][1:] <=55)
mask2 = mask2 #& mask3
avg_rollB = np.mean(r4B[1][1:][mask2])
# x_data2 = r4B[1][mask2]-avg_pitchB
x_data2 = r4B[1][1:][mask2]
y_data2 = -2.5*np.log(r4B[0][mask2])+2.5*np.log(np.max(r4B[0][mask2]))
# y_data2 = -2.5*np.log(r4B[0][mask2])
params_02B = curve_fit(harm,x_data2,y_data2,[130,18],maxfev=10000,bounds=((125,18),(350,23)))

fig_rfit,ax_rfit = plt.subplots(figsize=fig_size)
ax_rfit.plot(x_data,y_data,label=r'$U_{mem} (\phi)$',alpha=0.7,color='darkorange',linestyle='dashed',marker='o',markersize=m_size)
ax_rfit.plot(x_data,harm(x_data,*params_02A[0]),label=r'$U_{fit} (\phi)$',alpha=0.7,color='k')
ax_rfit.legend(fontsize=leg_size)
# ax_rfit.set_title(r'Roll ($\phi$) - Chain A')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02A[0][0])
# s2 = 'n = {:.2f}'.format(params_02A[0][1])
s2 = r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
# ax_rfit.text(np.mean(x_data),10,s1)
# ax_rfit.text(np.mean(x_data),12,s2)
# ax_rfit.text(0,14,s2)
fig_rfit.tight_layout()
ax_rfit.tick_params(labelsize=t_size)
ax_rfit.set_xlabel(r'$\phi$ (deg)',fontsize=f_size)
ax_rfit.set_ylabel("Potential Energy (kJ/mol)",fontsize=f_size)

fig_rfit,ax_rfit2 = plt.subplots(figsize=fig_size)
ax_rfit2.plot(x_data2,y_data2,label=r'$U_{mem} (\phi)$',alpha=0.7,color='darkorange',linestyle='dashed',marker='o',markersize=m_size)
ax_rfit2.plot(x_data2,harm(x_data2,*params_02B[0]),label=r'$U_{fit} (\phi)$',alpha=0.7,color='k')
ax_rfit2.legend(fontsize=leg_size)
# ax_rfit2.set_title(r'Roll ($\phi$) - Chain B')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02B[0][0])
# s2 = 'n = {:.2f}'.format(params_02B[0][1])
# ax_rfit2.text(np.mean(x_data2),10,s1)
s2 = r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitchB))
# ax_rfit2.text(np.mean(x_data2),12,s2)
# ax_rfit2.text(0,14,s2)
fig_rfit.tight_layout()
ax_rfit2.tick_params(labelsize=t_size)
ax_rfit2.set_xlabel(r'$\phi$ (deg)',fontsize=f_size)
ax_rfit2.set_ylabel("Potential Energy (kJ/mol)",fontsize=f_size)

print("Roll Data")
print("Avg for chainA : ",np.radians(avg_roll))
print("Avg for chainB: ", np.radians(avg_rollB))



k = np.linspace(1,100,10)
avg = [0.01,0.1,0.5,1,2]
x=np.radians(np.linspace(60,120,500))


fig,ax = plt.subplots()
for i in range(len(k)):
    u=0.5*k[i]*(1-np.cos(1*(x-90.0)))
    ax.plot(x,u,label=k[i],alpha=0.7,linewidth=1.0)

ax.legend()
ax.set_xlabel('r (nm)')
ax.set_ylabel('U(r)')


fig,ax = plt.subplots()

filtered_data = final_data[final_data['Timestep'].isin(filter_time)]
ax.hist(filtered_data['d1'],bins=100,label='d1',alpha=0.6,density=True,color='crimson')
ax.hist(filtered_data['d2'],bins=100,label='d2',alpha=0.6,density=True,color='royalblue')
ax.legend()   

plt.show()


fig,axA = plt.subplots(figsize=fig_size)

p4A=axA.hist(dataA['PitchA3'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='limegreen')
fig.tight_layout()
axA.tick_params(labelsize=t_size)
axA.set_xlabel(r'$\psi$ (deg)',fontsize=f_size)
axA.set_ylabel(r"$P_{mem} (\psi)$",fontsize=f_size)

fig,axB = plt.subplots(figsize=fig_size)
p4B=axB.hist(dataB['PitchB3'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='limegreen')
fig.tight_layout()
axB.tick_params(labelsize=t_size)
axB.set_xlabel(r'$\psi$ (deg)',fontsize=f_size)
axB.set_ylabel(r"$P_{mem} (\psi)$",fontsize=f_size)

fig,axA = plt.subplots(figsize=fig_size)
p4A=axA.hist(dataA['RollA4'][thresholdA],bins=100,label='chain A',alpha=0.6,density=True,color='darkorange')
fig.tight_layout()
axA.tick_params(labelsize=t_size)
axA.set_xlabel(r'$\phi$ (deg)',fontsize=f_size)
axA.set_ylabel(r"$P_{mem} (\phi)$",fontsize=f_size)

fig,axB = plt.subplots(figsize=fig_size)
p4B=axB.hist(dataB['RollB4'][thresholdB],bins=100,label='chain B',alpha=0.6,density=True,color='darkorange')
fig.tight_layout()
axB.tick_params(labelsize=t_size)
axB.set_xlabel(r'$\phi$ (deg)',fontsize=f_size)
axB.set_ylabel(r"$P_{mem} (\phi)$",fontsize=f_size)

plt.show()