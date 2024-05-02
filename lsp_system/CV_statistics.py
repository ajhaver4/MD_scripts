import pandas
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.optimize import curve_fit

#data = pandas.read_csv(sys.argv[1],delimiter="\t",comment='#',names=['frame','timestep','d1','d2'])

data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep='\s+',engine='python',names = ['time','d1','d2','m','n','q'])
#data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep='\s+',engine='python',names = ['time','d1','d2','m','n','q','r'])
#data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep=',',engine='python',names = ['time','d1','d2'])
# print(data['d1'][:10])
# sys.exit()
#Histogram

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

win_length = 400
dt = round(data['time'][1] - data['time'][0]) #ps
print("Time Step: ", dt, " ps")

#windows = np.arange(0,len(data['d1']),win_length)
#windows = (windows+win_length)*dt



d1_mean_vs_winsize = []
d1_std_vs_winsize = []
d2_mean_vs_winsize = []
d2_std_vs_winsize = []
windows = [5,10,20,50,100,200,400,500,1000,2000]

n = len(windows)
col = 3
rows = n/3

for win_length in windows:
    d1_means = []
    d1_std = []
    d2_means = []
    d2_std = []
    fig1,ax1 = plt.subplots()
    fig2,ax2 = plt.subplots()

    for i in range(0,len(data['d1']),win_length):
        hist_array = np.histogram(data['d1'],bins=50,density=True)
        print(hist_array[0])
        print(np.sum(hist_array[0]))

        window_time = (i+win_length)*dt
        plot_hist = ax1.hist(data['d1'][i:i+win_length],bins=50,density=True,alpha=0.5,label=str(window_time)+" ps",histtype='bar')
        bins_edges = plot_hist[1]
        bin_width = bins_edges[1]-bins_edges[0]
        print(bin_width)
        print(np.sum(plot_hist[0]*bin_width))

        #d1
        mean = np.mean(data['d1'][i:i+win_length])
        std = np.std(data['d1'][i:i+win_length])
        #param,pcov = curve_fit(gaus,plot_hist[1][:-1],plot_hist[0],p0=[1,mean,std])
        #plt.show()
        #sys.exit()


        plot_hist_d2 = ax2.hist(data['d2'][i:i+win_length],bins=50,density=True,alpha=0.5,label=str(window_time)+" ps",histtype='stepfilled')
        bins_edges_d2 = plot_hist_d2[1]
        bin_width_d2 = bins_edges_d2[1]-bins_edges_d2[0]
        print(bin_width_d2)
        print(np.sum(plot_hist_d2[0]*bin_width_d2))

        #d2
        mean_d2 = np.mean(data['d2'][i:i+win_length])
        std_d2 = np.std(data['d2'][i:i+win_length])
        #param_d2,pcov_d2 = curve_fit(gaus,plot_hist_d2[1][:-1],plot_hist_d2[0],p0=[1,mean_d2,std_d2])

        d1_means.append(mean)
        d1_std.append(std)

        d2_means.append(mean_d2)
        d2_std.append(std_d2)

    d1_mean_vs_winsize.append(np.mean(d1_means))
    d1_std_vs_winsize.append(np.mean(d1_std))

    d2_mean_vs_winsize.append(np.mean(d2_means))
    d2_std_vs_winsize.append(np.mean(d2_std))


    ax1_str = 'Distribution of d1 values w.r.t to ' + str(win_length*dt) + ' ps window size'
    ax2_str = 'Distribution of d2 values w.r.t to ' + str(win_length*dt) + ' ps window size'
    ax1.set_title(ax1_str)
    ax2.set_title(ax2_str)
    ax1.legend()
    ax2.legend()

    plt.show()




#d1
fig_3,ax_3 = plt.subplots()
plot_hist3 = ax_3.hist(data['d1'],bins=100,density=True)
mean = np.mean(data['d1'])
std = np.std(data['d1'])
param,pcov = curve_fit(gaus,plot_hist3[1][:-1],plot_hist3[0],p0=[1,mean,std])

#d2
fig_4,ax_4 = plt.subplots()
plot_hist4 = ax_4.hist(data['d2'],bins=100,density=True)
mean_d2 = np.mean(data['d2'])
std_d2 = np.std(data['d2'])
param_d2,pcov_d2 = curve_fit(gaus,plot_hist4[1][:-1],plot_hist4[0],p0=[1,mean_d2,std_d2])


ax_3.plot(plot_hist3[1][:-1],gaus(plot_hist3[1][:-1],*param),label='Gauss fitted d1')
ax3_str1 = r'$\mu = $' + str(round(mean,3))
ax3_str2 = r'$\sigma = $' + str(round(std,3))
ax_3.text(mean-std,max(plot_hist3[0]),ax3_str1)
ax_3.text(mean-std,max(plot_hist3[0])-1,ax3_str2)

ax_4.plot(plot_hist4[1][:-1],gaus(plot_hist4[1][:-1],*param_d2),label='Gauss fit d2')
ax4_str1 = r'$\mu = $' + str(round(mean_d2,3))
ax4_str2 = r'$\sigma = $' + str(round(std_d2,3))
ax_4.text(mean_d2-std_d2,max(plot_hist4[0]),ax4_str1)
ax_4.text(mean_d2-std_d2,max(plot_hist4[0])-1,ax4_str2)
ax_3.legend()
ax_4.legend()
#plt.show()


print("Gaussian fit parameters d1")
print(param)

print("Statistics of d1")
print("Mean: ")
print(np.mean(data['d1']))
print("STD: ")
print(np.std(data['d1']))

print("Gaussian fit parameters d2")
print(param_d2)

print("Statistics of d2")
print("Mean: ")
print(np.mean(data['d2']))
print("STD: ")
print(np.std(data['d2']))
#sys.exit()

# sig = np.std(data['d1'])
# sig_2 = np.std(data['d2'])
#
# n_frames = len(data['d1'])
# #delta_array = []
# print("Checking times for 1 sigma displacement")
#
#
# def check_times(data,sig,sig_2,n_frames):
#     d10 = data['d1'][0]
#     d20 = data['d2'][0]
#     time = 0
#     time_array = []
#
#     time_d2=0
#     time_array_d2 = []
#     delta_array = []
#
#     dt = 2 #ps
#
#     for i in range(1,n_frames):
#         dev = data['d1'][i] - d10
#         time=time+1
#
#         dev_2 = data['d2'][i] - d20
#         time_d2 = time_d2+1
#
#         if abs(dev) >= sig :
#             # print(i)
#             # print(data['d1'][i])
#             # print(d10)
#             d10 = data['d1'][i]
#             time_array.append(time*dt)
#             time=0
#
#         if abs(dev_2) >= sig_2 :
#             # print(i)
#             # print(data['d1'][i])
#             # print(d10)
#             d20 = data['d2'][i]
#             time_array_d2.append(time_d2*dt)
#             time_d2=0
#
#
#         delta_array.append(dev)
#     return(time_array,time_array_d2)
#
# (time_array,time_array_d2) = check_times(data,sig,sig_2,n_frames)
#
# print("d1 times")
# print(time_array)
# print(len(time_array))
# print("Average: ",np.mean(time_array[:-1]))
#
# print("d2 times")
# print(time_array_d2)
# print(len(time_array_d2))
# print("Average: ",np.mean(time_array_d2[:-1]))
#
# print("2 sigma")
# (time_array,time_array_d2) = check_times(data,0.1,0.1,n_frames)
# print("d1 times")
# print(time_array)
# print(len(time_array))
# print("Average: ",np.mean(time_array))
#
# print("d2 times")
# print(time_array_d2)
# print(len(time_array_d2))
# print("Average: ",np.mean(time_array_d2))

# fig,ax = plt.subplots()
#
# ax.plot(data['d1'],label='d1')
# ax.plot(data['d2'],label='d2')
# ax.legend()
# ax.set_xlabel("Frames")
# ax.set_ylabel("Distance (nm)")
window_times = np.array(windows)*dt
fig_mean,ax_mean = plt.subplots()
ax_mean.plot(window_times,np.array(d1_mean_vs_winsize),'o',label='Mean d1')
ax_mean.plot(window_times,np.array(d2_mean_vs_winsize),'o',label='Mean d2')

ax_mean.plot(window_times,np.array(d1_std_vs_winsize),'x', label = 'STD d1 ')
ax_mean.plot(window_times,np.array(d2_std_vs_winsize),'x',label = 'STD d2')
ax_mean.legend()
ax_mean.set_xlabel("Window length (ps)")
ax_mean.set_ylabel("nm")
ax_mean.set_title("Mean and std w.r.t window size")

plt.show()
print(d1_mean_vs_winsize)
print(d2_mean_vs_winsize)
