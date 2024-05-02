import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = pd.read_csv(sys.argv[1],comment='#',delimiter="\t")

n_frames = len(data['Timestep'])
timesteps_array = data['Timestep']/1e6
threshold = data['Timestep'] > 0e6
figp,[axp1,axp2,axp3,axp4,axp5] = plt.subplots(5,1)
p1A = axp1.hist(data['Z_angle_A'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p1B = axp1.hist(data['Z_angle_B'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp1.set_title("Z-angle")
axp1.legend()

p2A=axp2.hist(data['PitchA1'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p2B=axp2.hist(data['PitchB1'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp2.set_title("Pitch 1")
axp2.legend()
p3A=axp3.hist(data['PitchA2'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p3B=axp3.hist(data['PitchB2'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp3.set_title("Pitch 2")
axp3.legend()

p4A=axp4.hist(data['PitchA3'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
p4B=axp4.hist(data['PitchB3'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axp4.set_title("Pitch 3")
axp4.legend()

# p5A=axp5.hist(data['PitchA4'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
# p5B=axp5.hist(data['PitchB4'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
# axp5.set_title("Pitch 4")
# axp5.legend()

figr,[axr1,axr2,axr3,axr4,axr5] = plt.subplots(5,1)
r1A = axr1.hist(data['RollA1'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r1B = axr1.hist(data['RollB1'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr1.set_title("Roll 1")
axr1.legend()

r2A=axr2.hist(data['RollA2'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r2B=axr2.hist(data['RollB2'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr2.set_title("Roll 2")
axr2.legend()

r3A=axr3.hist(data['RollA3'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r3B=axr3.hist(data['RollB3'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr3.set_title("Roll 3")
axr3.legend()

r4A=axr4.hist(data['RollA4'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
r4B=axr4.hist(data['RollB4'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
axr4.set_title("Roll 4")
axr4.legend()

# r5A=axr5.hist(data['RollA5'][threshold],bins=100,label='chain A',alpha=0.6,density=True,color='crimson')
# r5B=axr5.hist(data['RollB5'][threshold],bins=100,label='chain B',alpha=0.6,density=True,color='gold')
# axr5.set_title("Roll 5")
# axr5.legend()


figA,[axA,axA1] = plt.subplots(2,1,sharex=True)
axA.plot(timesteps_array,data['Z_angle_A'],label='Z-angle',linewidth=0.7)
axA.plot(timesteps_array,data['PitchA1'],label='PitchA1',linewidth=0.7)
axA.plot(timesteps_array,data['PitchA2'],label='PitchA2',linewidth=0.7)
axA.plot(timesteps_array,data['PitchA3'],label='PitchA3',linewidth=0.7)
# axA.plot(timesteps_array,data['PitchA4'],label='PitchA4',linewidth=0.7)
axA1.plot(timesteps_array,data['Z_angle_B'],label='Z-angle',linewidth=0.7)
axA1.plot(timesteps_array,data['PitchB1'],label='PitchB1',linewidth=0.7)
axA1.plot(timesteps_array,data['PitchB2'],label='PitchB2',linewidth=0.7)
axA1.plot(timesteps_array,data['PitchB3'],label='PitchB3',linewidth=0.7)
# axA1.plot(timesteps_array,data['PitchB4'],label='PitchB4',linewidth=0.7)
axA.legend()
axA1.legend()
axA.set_xlabel("Time")
axA.set_ylabel("Angle ")
axA1.set_ylabel("Angle ")
axA.set_title("Chain A")
axA1.set_title("Chain B")

figB,[axB,axB1] = plt.subplots(2,1,sharex=True)
axB.plot(timesteps_array,data['RollA1'],label="Roll1",linewidth=0.7)
axB.plot(timesteps_array,data['RollA2'],label="Roll2",linewidth=0.7)
axB.plot(timesteps_array,data['RollA3'],label="Roll3",linewidth=0.7)
axB.plot(timesteps_array,data['RollA4'],label="Roll4",linewidth=0.7)
# axB.plot(timesteps_array,data['RollA5'],label="Roll5",linewidth=0.7)
axB1.plot(timesteps_array,data['RollB1'],label="Roll1",linewidth=0.7)
axB1.plot(timesteps_array,data['RollB2'],label="Roll2",linewidth=0.7)
axB1.plot(timesteps_array,data['RollB3'],label="Roll3",linewidth=0.7)
axB1.plot(timesteps_array,data['RollB4'],label="Roll4",linewidth=0.7)
# axB1.plot(timesteps_array,data['RollB5'],label="Roll5",linewidth=0.7)
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
    v = 0.5*k*(np.radians(theta-theta0))**2 #- 2.5*np.log(c)
    return(v)

"""Chain A Pitch"""

hist_pitchA = p4A

mask1 = np.nonzero(hist_pitchA[0])
avg_pitch = np.mean(p4A[1][mask1])
# x_data = p4A[1][mask1]-avg_pitch
x_data = hist_pitchA[1][mask1]
y_data = -2.5*np.log(hist_pitchA[0][mask1])+2.5*np.log(np.max(hist_pitchA[0][mask1]))
# y_data = -2.5*np.log(hist_pitchA[0][mask1])
params_02 = curve_fit(harm,x_data,y_data,[100,90],maxfev=8000,bounds=(1.0,2000))
# params_02 = curve_fit(harm2,x_data,y_data,400)
#
#
#
fig_fit,[ax_fit,ax_fit2] = plt.subplots(2,1)
ax_fit.plot(x_data,y_data,label='V(x)=-2.5 * log(P(x))')
ax_fit.plot(x_data,harm(x_data,*params_02[0]),label='V(x) = k*(1-cos(theta))')
ax_fit.legend()
ax_fit.set_title('PitchA3')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02[0][0])
# s2 = 'n = {:.2f}'.format(params_02[0][1])
s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
ax_fit.text(np.mean(x_data),10,s1)
ax_fit.text(np.mean(x_data),12,s2)
# ax_fit.text(0,14,s2)

"""Chain B Pitch"""

hist_pitchB = p4B
mask2 = np.nonzero(hist_pitchB[0])
avg_pitchB = np.mean(p4B[1][mask2])
# x_data1 = p4B[1][mask2]-avg_pitchB
x_data1 = hist_pitchB[1][mask2]
y_data1 = -2.5*np.log(hist_pitchB[0][mask2]) + 2.5*np.log(np.max(hist_pitchB[0][mask2]))
# y_data1 = -2.5*np.log(hist_pitchB[0][mask2])
print("Pitch Shift U(mean): ", 2.5*np.log(np.max(hist_pitchB[0][mask2])))
params_02B = curve_fit(harm,x_data1,y_data1,[100,90],maxfev=8000,bounds=(1.0,2000))
# params_02B = curve_fit(harm2,x_data1,y_data1,400,maxfev=800,bounds=(1.0,2000))

print("Pitch Data")
print("Avg for chainA : ",np.radians(avg_pitch))
print("Avg for chainB: ", np.radians(avg_pitchB))


ax_fit2.plot(x_data1,y_data1,label='V(x)=-2.5 * log(P(x))')
# ax_fit2.plot(x_data1,harm(x_data1,params_02[0][0]),label='V(x) = k*(1-cos(theta)) for A')
ax_fit2.plot(x_data1,harm(x_data1,*params_02B[0]),label='V(x) = k*(1-cos(theta))')
ax_fit2.legend()
ax_fit2.set_title('PitchB3')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02B[0][0])
# s2 = 'n = {:f}'.format(params_02B[0][1])
s2=r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitchB))
ax_fit2.text(np.mean(x_data1),10,s1)
ax_fit2.text(np.mean(x_data1),12,s2)
# ax_fit2.text(0,14,s2)

"""Chain A Roll"""
mask1 = np.nonzero(r4A[0])
avg_pitch = np.mean(data['RollA4'][threshold])
# x_data = r4A[1][mask1]-avg_pitch
x_data = r4A[1][mask1]
y_data = -2.5*np.log(r4A[0][mask1])+2.5*np.log(np.max(r4A[0][mask1]))
# y_data = -2.5*np.log(r4A[0][mask1])
params_02A = curve_fit(harm,x_data,y_data,[230,18])

"""Chain B Roll"""
mask2 = np.nonzero(r4B[0])
avg_pitchB = np.mean(data['RollB4'][threshold])
# x_data2 = r4B[1][mask2]-avg_pitchB
x_data2 = r4B[1][mask2]
y_data2 = -2.5*np.log(r4B[0][mask2])+2.5*np.log(np.max(r4B[0][mask2]))
# y_data2 = -2.5*np.log(r4B[0][mask2])
params_02B = curve_fit(harm,x_data2,y_data2,[240,25],maxfev=10000)

fig_rfit,[ax_rfit,ax_rfit2] = plt.subplots(2,1)
ax_rfit.plot(x_data,y_data,label='V(x)=-2.5 * log(P(x))')
ax_rfit.plot(x_data,harm(x_data,*params_02A[0]),label='V(x) = k*(1-cos(theta))')
ax_rfit.legend()
ax_rfit.set_title('RollA4')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02A[0][0])
# s2 = 'n = {:.2f}'.format(params_02A[0][1])
s2 = r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitch))
ax_rfit.text(np.mean(x_data),10,s1)
ax_rfit.text(np.mean(x_data),12,s2)
# ax_rfit.text(0,14,s2)

ax_rfit2.plot(x_data2,y_data2,label='V(x)=-2.5 * log(P(x))')
ax_rfit2.plot(x_data2,harm(x_data2,*params_02B[0]),label='V(x) = k*(1-cos(theta))')
ax_rfit2.legend()
ax_rfit2.set_title('RollB4')
s1 = 'k = {:.2f} kJ mol-1 rad-2'.format(params_02B[0][0])
# s2 = 'n = {:.2f}'.format(params_02B[0][1])
ax_rfit2.text(np.mean(x_data2),10,s1)
s2 = r'$\mu$ = {:.2f} rad'.format(np.radians(avg_pitchB))
ax_rfit2.text(np.mean(x_data2),12,s2)
# ax_rfit2.text(0,14,s2)

print("Roll Data")
print("Avg for chainA : ",np.radians(avg_pitch))
print("Avg for chainB: ", np.radians(avg_pitchB))




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
plt.show()
