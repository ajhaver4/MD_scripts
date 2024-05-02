import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

main_df = pd.DataFrame()
rmsf_A = None
rmsd_B = None
fig,[ax1,ax2] = plt.subplots(2,1)
for file in sys.argv[1:]:
    time=file.split("_")[-1]

    data = pd.read_csv(file,comment='#',delimiter="\t")
    if rmsf_A is None:
        rmsf_A = np.reshape(data['RMSFA'].values,(len(data['ResidueA']),1))
        rmsf_B = np.reshape(data['RMSFB'].values,(len(data['ResidueA']),1))
    else:

        rmsf_A = np.concatenate((rmsf_A,np.reshape(data['RMSFA'].values,(len(data['ResidueA']),1))),axis=1)
        rmsf_B = np.hstack((rmsf_B,np.reshape(data['RMSFB'].values,(len(data['ResidueA']),1))))


    ax1.plot(data["ResidueA"],data["RMSFA"],linewidth=0.8,label=time)
    ax2.plot(data["ResidueB"],data['RMSFB'],linewidth=0.8,label=time)
    main_df=pd.concat([main_df,data])



avg_rmsfA = np.mean(rmsf_A,axis=1)
avg_rmsfB = np.mean(rmsf_B,axis=1)
ax1.plot(data["ResidueA"],avg_rmsfA,linewidth=1.0,linestyle='dashed',label="mean",color='black')
ax2.plot(data["ResidueB"],avg_rmsfB,linewidth=1.0,linestyle='dashed',label="mean",color='black')

ax1.set_title("Chain A RMSF")
ax2.set_title("Chain B RMSF")
ax1.legend()
ax2.legend()
plt.show()

#Identify regions of high flexibility
high_flexible_A = []
low_flexible_A = []


maskA_1 = avg_rmsfA <= 1.0
maskA_2 = avg_rmsfA >5

maskB_1 = avg_rmsfB <=1.0
maskB_2 = avg_rmsfB >5

print("------- Chain A ------")
print("Atoms with high fluctuations: ", list(data['ResidueA'][maskA_2].values))
print("Atoms with low fluctuations: ",list(data['ResidueA'][maskA_1].values))

print("----- Chain B--------")

print("Atoms with high fluctuations: ", list(data['ResidueB'][maskB_2].values))
print("Atoms with low fluctuations: ",list(data['ResidueB'][maskB_1].values))
