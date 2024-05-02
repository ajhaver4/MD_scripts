import sys
import numpy as np
import pandas
import matplotlib.pyplot as plt

#data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep='\s+',engine='python',names = ['time','d1','d2'])
#data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep='\s+',engine='python',names = ['time','d1','d2','bias','rct','rbias','dist1'])
data = pandas.read_csv(sys.argv[1],comment='#',header=None,sep=',',engine='python',names = ['time','d1','d2'])


mean_d1 = []
mean_d2 = []
time_lag = []
lag = 10
while lag < len(data['time']) -1 :
    avg_delta_d1 = []
    avg_delta_d2 = []
    iter=0
    print("Lag time: ",lag)
    while iter < len(data['time']) - lag:
        delta_d1 = data['d1'][iter+lag] - data['d1'][iter]
        delta_d2 = data['d2'][iter+lag] - data['d2'][iter]

        avg_delta_d1.append(abs(delta_d1))
        avg_delta_d2.append(abs(delta_d2))

        iter+=1
    print(len(avg_delta_d1))
    mean_d1.append(np.mean(avg_delta_d1))
    mean_d2.append(np.mean(avg_delta_d2))

    time_lag.append(lag)
    lag+=1

#Second method
# while lag < len(data['time']) -1 :
#     avg_delta_d1 = []
#     avg_delta_d2 = []
#     iter=0
#     print("Lag time: ",lag)
#     while iter < len(data['time']) - lag:
#         delta_d1 = data['d1'][iter+lag] - data['d1'][iter]
#         delta_d2 = data['d2'][iter+lag] - data['d2'][iter]
#         avg_delta_d1.append(abs(delta_d1))
#         avg_delta_d2.append(abs(delta_d2))
#
#         iter+=lag
#     print(len(avg_delta_d1))
#     mean_d1.append(np.mean(avg_delta_d1))
#     mean_d2.append(np.mean(avg_delta_d2))
#
#     time_lag.append(lag)
#     lag+=1


fig,ax = plt.subplots()
ax.plot(time_lag,mean_d1,label='d1',linewidth=0.6)
ax.plot(time_lag,mean_d2,label='d2',linewidth=0.6)
ax.legend()
ax.set_xlabel("Frames ")
ax.set_ylabel("Mean change in CVs")
plt.show()
