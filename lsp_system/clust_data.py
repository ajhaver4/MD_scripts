import pandas
import sys
import matplotlib.pyplot as plt
import math
import numpy as np
from sklearn.cluster import KMeans
from matplotlib.markers import MarkerStyle as markerstyle
from mpl_toolkits.mplot3d import axes3d

file=sys.argv[1]
num_clusters = int(sys.argv[2])
cut_off = float(sys.argv[3])
border = float(sys.argv[4])
data = pandas.read_csv(file,comment='#',sep='\t')

n_values = len(data['d1'])
#print(data['Timestep'][3270:3280])

"""
Clustering
"""
def cluster_method(n_cluster,dist_mat):
    #KMEANS Clustering
    kmeans = KMeans(n_clusters = n_cluster, random_state=0).fit(dist_mat)
    clus_cen = kmeans.cluster_centers_ #Obtain centroids for all the clusters
    transform_mat = kmeans.transform(dist_mat) #This method calculates the distance of each point from each cluster
    labels = kmeans.labels_  #Labels which frame belongs to which cluster
    #print("KMeans labels: ",kmeans.labels_)
    #print("Cluster centers:" ,kmeans.cluster_centers_)
    return(labels,clus_cen,transform_mat)

dist_mat2 = []
rmsd_mat = []
angle_dist_mat = []
langle_dist_mat = []
min_dist = []
timesteps = []
count = 0.0
print("Selecting frames as per cut-off...")
for i in range(n_values):
    if data['Min dist'][i] < cut_off and data['d1'][i] + data['d2'][i] < border:
        mean_d = (data['d1'][i]+data['d2'][i])/2.0
        dist_mat2.append([data['d1'][i],data['d2'][i]])
        rmsd_mat.append([mean_d,data['RMSD-B'][i]])
        #angle_dist_mat.append(mean_d[i],data['ang'][i]])
        langle_dist_mat.append([mean_d,data['Ax-ang'][i]])
        min_dist.append(data['Min dist'][i])
        timesteps.append(data['Timestep'][i])
        count+=1

print("Number of selected frames: ",count)
print("Clustering performed on {:f} frames of the trajectory.".format(count*100/n_values))
dist_mat2 = np.reshape(np.array(dist_mat2),(len(dist_mat2),2))
rmsd_mat = np.reshape(np.array(rmsd_mat),(len(rmsd_mat),2))
#angle_dist_mat = np.reshape(np.array(angle_dist_mat),(len(angle_dist_mat),2))
langle_dist_mat = np.reshape(np.array(langle_dist_mat),(len(langle_dist_mat),2))

#Dist_CV Clustering
print("Begin distance clustering...")
(labels,clus_cen,transform_dist_mat) = cluster_method(num_clusters,dist_mat2)
clustal_data = {ind : {} for ind in range(num_clusters)}
timestep_data = {ind : [] for ind in range(num_clusters)}
for i in range(len(labels)):
    cluster_index = labels[i]
    dist_to_clus = transform_dist_mat[i][cluster_index]
    clustal_data[cluster_index][dist_to_clus] = i
    timestep_data[cluster_index].append(timesteps[i])

#Need to find the closest frame to a cluster centorid
#This frame would be the reference to calculate the RMSD
#This frame would then be used for the ensemble image presentation
min_dist_frames = {i:0 for i in range(num_clusters)}
for i in range(num_clusters):
    distances = list(clustal_data[i].keys())
    #print(distances)
    keys = clustal_data.keys()
    #print("Keys for clustal_data: ",type(keys))
    min_dist_frames[list(clustal_data.keys())[i]] = clustal_data[i][(min(distances))]
"""
#Perform clustering based on RMSD and Min_distance:
print("RMSD Clustering...")
(rmsd_labels,rmsd_clus_cen,rmsd_transform_dist_mat) = cluster_method(num_clusters,rmsd_mat)
rmsd_clustal_data ={ind : {} for ind in range(num_clusters)}
for i in range(len(rmsd_labels)):
    cluster_index = rmsd_labels[i]
    dist_to_clus = rmsd_transform_dist_mat[i][cluster_index]
    rmsd_clustal_data[cluster_index][dist_to_clus] = i

print("Begin Angle Clustering")
#Clustering based on angle and dist CVs

(angle_labels,angle_clus_cen,transform_angle_mat) = cluster_method(num_clusters,angle_dist_mat)
angle_clustal_data = {ind : {} for ind in range(num_clusters)}
for i in range(len(angle_labels)):
    cluster_index = angle_labels[i]
    dist_to_clus = transform_angle_mat[i][cluster_index]
    angle_clustal_data[cluster_index][dist_to_clus] = i

#Clustering based on axial angle and dist
(langle_labels,langle_clus_cen,transform_langle_mat) = cluster_method(num_clusters,langle_dist_mat)
langle_clustal_data = {ind : {} for ind in range(num_clusters)}
for i in range(len(langle_labels)):
    cluster_index = langle_labels[i]
    dist_to_clus = transform_langle_mat[i][cluster_index]
    langle_clustal_data[cluster_index][dist_to_clus] = i
"""
"""
#--------------------------------------------------------------------------------
#--------------------------- PLOTTING/WRITING TRAJ -------------------------------------------
#--------------------------------------------------------------------------------
"""
with open("Cluster_results_rmsd.txt", "w") as fl:
    fl.write("Clustering performed with :")
    fl.write("Cut-off: ")
    fl.write(str(sys.argv[3]))
    fl.write("\tBorder: ")
    fl.write(str(sys.argv[4]))
    fl.write("\n")
    fl.write("#Cluster centroids")
    fl.write("\n")
    fl.write("#Cluster No.")
    fl.write("\t")
    fl.write("x")
    fl.write("\t")
    fl.write("y")
    fl.write("\n")
    for i in range(len(clus_cen)):
        fl.write(str(i))
        fl.write("\t")
        fl.write(str(clus_cen[i,0]))
        fl.write("\t")
        fl.write(str(clus_cen[i,1]))
        fl.write("\n")
        d1_frame = np.where(data['d1'] == dist_mat2[min_dist_frames[i],0])
        d2_frame = np.where(data['d2'] == dist_mat2[min_dist_frames[i],1])
        fl.write("Minimum distance frame: ")
        fl.write("\t")
        fl.write("Timestep: ")
        fl.write("\t")
        fl.write(str(data['Timestep'][d1_frame[0]]))
        fl.write("\n")
        fl.write("Minimum distance frame CV Coords: ")
        fl.write("\t")
        fl.write(str(dist_mat2[min_dist_frames[i],0]))
        fl.write("\t")
        fl.write(str(dist_mat2[min_dist_frames[i],1]))
        fl.write("\n")
        fl.write("Dist b/w memb bind residues: ")
        fl.write(str(data['Memb dist'][d1_frame[0]]))
        fl.write("\n")
        # fl.write("Rot angle: ")
        # fl.write(str(data['Rot_angle'][d1_frame[0]]))
        # fl.write("\n")
        fl.write("RMSD: ")
        fl.write(str(data['RMSD-B'][d1_frame[0]]))
        fl.write("\n")
        fl.write("Number of frames: ")
        fl.write(str(len(list(clustal_data[i].values()))))
        fl.write("\n")
        fl.write("Frame list: \n")
        fl.write(str(list(clustal_data[i].values())))
        fl.write("\n\n")

with open("Cluster_frames.txt", "w") as fl1:
    fl1.write("Clustering performed with :")
    fl1.write("Cut-off: ")
    fl1.write(str(sys.argv[3]))
    fl1.write("\tBorder: ")
    fl1.write(str(sys.argv[4]))
    fl1.write("\n")
    fl1.write("#Cluster No.")
    fl1.write("\n")
    for i in range(len(clus_cen)):
        fl1.write("@")
        fl1.write(str(i))
        fl1.write("\t")
        fl1.write("Timesteps: \n")
        fl1.write(str(timestep_data[i]))
        fl1.write("\n\n")

#Plotting the cluster representation
#Each frame is represented by a coloured coded circle depending upon the cluster it belongs to
#Cluster centroid are represented by '+'
#Closest frame to each cluster in represented by ^
#label_colours = {i:float(i) for i in range(num_clusters)}
label_colours = {0:'forestgreen', 1:'teal', 2:'crimson', 3:'gold', 4:'orchid', 5:'peru', 6:'mediumpurple', 7:'darkorange'}
sub_label_colours = {0:'crimson',1:'olivedrab',2:'teal',3:'darkorange'}
#c_map = [label_colours[i] for i in labels]
#print(labels.astype(np.float))
with open("centroids_"+str(num_clusters)+"_04.txt", 'w') as fl3:
    fl3.write("#Defining cluster centroid as different bound states\n")
    fl3.write("States\t")
    fl3.write("xmin\txmax\tymin\tymax\n")
    for i in range(len(clus_cen)):
        fl3.write("C")
        fl3.write(str(i))
        fl3.write("\t")
        fl3.write(str(dist_mat2[min_dist_frames[i],0]))
        fl3.write("\t")
        fl3.write(str(dist_mat2[min_dist_frames[i],0]))
        fl3.write("\t")
        fl3.write(str(dist_mat2[min_dist_frames[i],1]))
        fl3.write("\t")
        fl3.write(str(dist_mat2[min_dist_frames[i],1]))
        fl3.write("\n")
#
fig,ax = plt.subplots()
min_dist_frame_cc = []
for i in range(num_clusters):
    fr_num = list(clustal_data[i].values())
    plot_dist = np.zeros((len(fr_num),2))
    for j in range(len(fr_num)):
        plot_dist[j,:] = dist_mat2[fr_num[j],:2]
        if fr_num[j] == min_dist_frames[i]:
            min_dist_frame_cc.append(j)
    #c_map = np.ones(len(fr_num))*i
    #print(c_map)
    ax.scatter(plot_dist[:,0],plot_dist[:,1],c=label_colours[i],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
    #ax.scatter(plot_dist[:,0],plot_dist[:,1],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
ax.scatter(clus_cen[:,0],clus_cen[:,1],c='black',marker='+',label ='Centroids')
ax.legend(fancybox=True,framealpha=0.4,fontsize='large',handletextpad=0.1,markerscale=2.5,loc='best',ncol=num_clusters+1,columnspacing=0.5)
ax.set_title("d1 d2 Clustering")
ax.set_xlabel('d1 (nm)',fontsize='x-large')
ax.set_ylabel('d2 (nm)',fontsize='x-large')

fig.savefig("Bound Distance Clustering - "+str(num_clusters) + "_" + str(cut_off)+ " Clusters.png")

#plt.show()
"""
#Create another plot for sub cluster
fig5,ax5 = plt.subplots()
min_dist_frame_cc = []
for i in range(num_clusters):
    fr_num = list(clustal_data[i].values())
    plot_dist = np.zeros((len(fr_num),2))
    for j in range(len(fr_num)):
        plot_dist[j,:] = dist_mat2[fr_num[j],:2]
        if fr_num[j] == min_dist_frames[i]:
            min_dist_frame_cc.append(j)
    #c_map = np.ones(len(fr_num))*i
    #print(c_map)
    #ax5.scatter(plot_dist[:,0],plot_dist[:,1],c=label_colours[i],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
    ax5.scatter(plot_dist[:,0],plot_dist[:,1],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
ax5.scatter(clus_cen[:,0],clus_cen[:,1],c='black',marker='+',label ='Centroids')


for i in range(sub_num_clusters):
    fr_num = list(sub_clustal_data[i].values())
    sub_plot_dist = np.zeros((len(fr_num),2))
    for j in range(len(fr_num)):
        sub_plot_dist[j,:] = dist_mat2[fr_num[j],:]
    ax5.scatter(sub_plot_dist[:,0],sub_plot_dist[:,1],c=sub_label_colours[i],marker='*',label='Sub-Cluster: '+str(i))

ax5.scatter(sub_clus_cen[:,0],sub_clus_cen[:,1],c='black',marker='*',label='Sub-centroids')

ax5.legend(fancybox=True,framealpha=0.4,fontsize='large',handletextpad=0.1,markerscale=2.5,loc='best',ncol=num_clusters+1,columnspacing=0.5)
ax5.set_title("Sub Cluster")
ax5.set_xlabel('d1 (nm)',fontsize='x-large')
ax5.set_ylabel('d2 (nm)',fontsize='x-large')

fig5.savefig("Sub Clustering - "+str(num_clusters)+" Clusters.png")
"""

"""
fig3,ax3 = plt.subplots()
for i in range(num_clusters):
    fr_num = list(rmsd_clustal_data[i].values())
    plot_dist = np.zeros((len(fr_num),2))
    for j in range(len(fr_num)):
        plot_dist[j,:] = rmsd_mat[fr_num[j],:]

    #ax3.scatter(plot_dist[:,0],plot_dist[:,1],c=label_colours[i],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
    ax3.scatter(plot_dist[:,0],plot_dist[:,1],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
#ax3.legend(fancybox=True,framealpha=0.4,fontsize='large',handletextpad=0.1,markerscale=2.5,loc='best',ncol=num_clusters+1,columnspacing=0.5)
ax3.set_title("KMeans Cluster Analysis: RMSD")
ax3.set_xlabel('RMSD (nm)',fontsize='x-large')
ax3.set_ylabel('dist (nm)',fontsize='x-large')

fig3.savefig("RMSD Distance Clustering - "+str(num_clusters)+" Clusters.png")

fig4,ax4 = plt.subplots()
for i in range(num_clusters):
    fr_num_rms = list(rmsd_clustal_data[i].values())
    #fr_num = list(clustal_data[i].values())
    plot_dist = np.zeros((len(fr_num_rms),2))
    for j in range(len(fr_num_rms)):
        plot_dist[j,:] = dist_mat2[fr_num_rms[j],:2]
    ax4.scatter(plot_dist[:,0],plot_dist[:,1],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
    #ax4.scatter(plot_dist[:,0],plot_dist[:,1],c=label_colours[i],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
#ax4.legend(fancybox=True,framealpha=0.4,fontsize='large',handletextpad=0.1,markerscale=2.5,loc='best',ncol=num_clusters+1,columnspacing=0.5)
ax4.set_title("plot d1 d2 : cluster: RMSD")
ax4.set_xlabel('d1 (nm)',fontsize='x-large')
ax4.set_ylabel('d2 (nm)',fontsize='x-large')
fig4.savefig("RMSD Clustering d1 d2 plotting- "+str(num_clusters)+" Clusters.png")


#Plotting dist and ang cluster
fig_ang,ax_ang = plt.subplots()
for i in range(num_clusters):
    fr_num = list(angle_clustal_data[i].values())
    plot_angle = np.zeros((len(fr_num),2))
    for j in range(len(fr_num)):
        plot_angle[j,:] = angle_dist_mat[fr_num[j],:]
    #ax_ang.scatter(plot_angle[:,0],plot_angle[:,1],c=label_colours[i],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
    ax_ang.scatter(plot_angle[:,0],plot_angle[:,1],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
ax_ang.scatter(angle_clus_cen[:,0],angle_clus_cen[:,1],c='black',marker='+',label ='Centroids')
#ax_ang.legend(fancybox=True,framealpha=0.4)
ax_ang.set_title("Angle Clustering ang vs dist")
ax_ang.set_xlabel("dist (nm))")
ax_ang.set_ylabel("Angle")
fig_ang.savefig("Angle Clustering_01 - "+str(num_clusters)+" Clusters.png")


#Axial ANgle Clustering
fig_lang,ax_lang = plt.subplots()
for i in range(num_clusters):
    fr_num = list(langle_clustal_data[i].values())
    plot_langle = np.zeros((len(fr_num),2))
    for j in range(len(fr_num)):
        plot_langle[j,:] = langle_dist_mat[fr_num[j],:]
    #ax_lang.scatter(plot_langle[:,0],plot_langle[:,1],c=label_colours[i],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
    ax_lang.scatter(plot_langle[:,0],plot_langle[:,1],marker=markerstyle(marker='.',fillstyle='none'),label= 'C-'+str(i))
ax_lang.scatter(langle_clus_cen[:,0],langle_clus_cen[:,1],c='black',marker='+',label ='Centroids')
#ax_lang.legend(fancybox=True,framealpha=0.4)
ax_lang.set_title("Angle Clustering ang vs dist")
ax_lang.set_xlabel("dist (nm))")
ax_lang.set_ylabel("Angle")
fig_lang.savefig("Axial Angle Clustering_01 - "+str(num_clusters)+" Clusters.png")

fig_min,ax_min = plt.subplots()
ax_min.plot(min_dist,linewidth=0.6)
ax_min.set_ylabel('Distance nm')
ax_min.set_xlabel('Frame')
plt.savefig("Min_dist.png",dpi=300,bbox_inches='tight')
"""
plt.show()

#Plotting all variables in each cluster


for i in range(num_clusters):
    fr_num = list(clustal_data[i].values())
    rmsd_val = []
    dim_angle = []
    axial_ang = []
    memb_dist = []
    d3 = []
    fr_num.sort()
    for j in range(len(fr_num)):
        rmsd_val.append(data['RMSD-B'][fr_num[j]])
        #dim_angle.append(data['Dim_ang-1'][fr_num[j]])
        axial_ang.append(math.degrees(data['Ax-ang'][fr_num[j]]))
        #d3.append(data['dist'][fr_num[j]])

    fig_01,ax_01 = plt.subplots(3,1,sharex=True)
    ax_01[0].plot(d3,'r',linewidth=0.6)
    ax_01[1].plot(rmsd_val,'y',linewidth=0.6)
    ax_01[2].plot(axial_ang,'b',linewidth=0.6,label='Ax-ang')
    ax_01[2].plot(dim_angle,'forestgreen',linewidth=0.6,label='dim angle')
    ax_01[2].set_xlabel('Frame')
    ax_01[0].set_ylabel('nm')
    ax_01[1].set_ylabel('nm')
    ax_01[2].set_ylabel('degrees')
    ax_01[0].set_title('d3')
    ax_01[1].set_title('RMSD')
    ax_01[2].set_title('Angles')
    ax_01[2].legend()
    fig_01.suptitle('Cluster - ' + str(i))

    plt.show()

fig_dis,ax_dis = plt.subplots(4,1)
for i in range(num_clusters):
    fr_num = list(clustal_data[i].values())
    rmsd_val = []
    #ins_angle = []
    axial_ang = []
    memb_dist = []
    #d3 = []
    fr_num.sort()
    for j in range(len(fr_num)):
        rmsd_val.append(data['RMSD-B'][fr_num[j]])
        #ins_angle.append(data['ins-angle'][fr_num[j]])
        axial_ang.append(data['Ax-ang'][fr_num[j]])
        memb_dist.append(data['Memb dist'][fr_num[j]])


    #ax_01[0].plot(d3,'r',linewidth=0.6)
    lb_str = "C-" + str(i)
    ax_dis[0].hist(rmsd_val,bins=100,density=True,histtype='stepfilled',color=label_colours[i],alpha=0.6,label=lb_str)
    ax_dis[1].hist(axial_ang,bins=100,density=True,histtype='stepfilled',color=label_colours[i],alpha=0.6,label=lb_str)
    #ax_dis[2].hist(ins_angle,bins=100,density=True,histtype='stepfilled',color=label_colours[i],alpha=0.4,label=lb_str)
    ax_dis[2].hist(memb_dist,bins=100,density=True,histtype='stepfilled',color=label_colours[i],alpha=0.6,label=lb_str)
    #ax_01[2].plot(dim_angle,'forestgreen',linewidth=0.6,label='dim angle')

    #ax_dis[1].legend()
    #ax_dis[2].legend()
    #ax_dis[3].legend()
    #fig_01.suptitle('Cluster - ' + str(i))


ax_dis[1].axvline(x=142.74,linestyle='--',linewidth=0.8)
ax_dis[2].axvline(x=3.356,linestyle='--',linewidth=0.8)
#ax_dis[3].axvline(x=0.9328,linestyle='--',linewidth=0.8)
ax_dis[0].set_ylabel('Frequency')
ax_dis[1].set_ylabel('Frequency')
ax_dis[2].set_ylabel('Frequency')
#ax_dis[3].set_ylabel('Frequency')
ax_dis[0].set_title('RMSD')
ax_dis[1].set_title('Long Angle')
#ax_dis[2].set_title('Ins Helix Angle')
ax_dis[2].set_title('Memb distance')
ax_dis[0].legend(fancybox=True,framealpha=0.4,fontsize='medium',handletextpad=0.1,markerscale=2.5,loc='best',ncol=num_clusters+1,columnspacing=0.5)
fig_dis.tight_layout()
plt.show()


"""
print("Minimum Distance frames")
print(len(data['Timestep']))
for i in range(len(min_dist_frames)):
    print(min_dist_frames[i])
    print("Dist coords from data")
    print(data['d1'][min_dist_frames[i]])
    print(data['d2'][min_dist_frames[i]])
    print("dist_mat2 values")
    print(dist_mat2[min_dist_frames[i],0])
    print(dist_mat2[min_dist_frames[i],1])
    print("Timestep")
    d1_frame = np.where(data['d1'] == dist_mat2[min_dist_frames[i],0])
    d2_frame = np.where(data['d2'] == dist_mat2[min_dist_frames[i],1])
    print(data['Timestep'][d1_frame[0]])
    print("Frame: ",d1_frame[0])
"""
