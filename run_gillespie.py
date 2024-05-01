import numpy as np
import matplotlib.pyplot as plt
import random
import math
import sys

#Define starting parameter

#Molecules
n_A = 0
n_B = 0

#Rate constants
k1 = 1e-3           # A + A -> C
k2 = 1e-2           # A + B -> D
k3 = 1.2            # null -> A
k4 = 1             # null -> B

#Reaction Matrix
M = np.array([[-2,0,1,0],[-1,-1,0,1],[1,0,0,0],[0,1,0,0]])
timeConc = np.reshape(np.array([0,0,0,0]),(1,4))
currConc = np.reshape(np.array([0,0,0,0]),(1,4))
propMat = np.array([0.0,0.0,0.0,0.0])
reactIdx = 0
countRxn = 0

print(currConc.shape)

#Total time
totalTime = 100000  #sec
time = 0.0
tau = 0.0
time_array = [0.0]

numRxn = M.shape[0]


while time < totalTime:
    #Generate two random numbers
    r1 = np.random.uniform()
    r2 = np.random.uniform()

    #Calculate propensities
    propMat[0] = currConc[0,0]*(currConc[0,0]-1)*k1
    #print(propMat)
    propMat[1] = currConc[0,0]*currConc[0,1]*k2
    propMat[2] = k3
    propMat[3] = k4

    sum_prop = np.sum(propMat)

    #Calculate time of next reaction
    tau = (1/sum_prop)*math.log(1/r1)

    if time+tau > totalTime:
        break

    #print("Time of next reaction :",tau)

    #Calculate which reaction occurs
    cum_sum = 0
    for j in range(len(propMat)):
        cum_sum = cum_sum+(propMat[j]/sum_prop)

        if r2 < cum_sum:
            reactIdx = j
            #print("Reaction: ",j)
            break

    currConc = currConc + M[reactIdx,:]
    timeConc = np.concatenate((timeConc,currConc),axis=0)
    time = time+tau
    time_array.append(time)
    countRxn+=1

print("Number of reactions:", countRxn)
print(currConc)
#print(timeConc)

with open("conc_vs_time",'a') as fl:
    fl.write("#Time\tA\tB\n")
    for v in range(len(time_array)):
        fl.write(str(time_array[v])+"\t")
        fl.write(str(timeConc[v,0])+"\t")
        fl.write(str(timeConc[v,1])+"\n")
