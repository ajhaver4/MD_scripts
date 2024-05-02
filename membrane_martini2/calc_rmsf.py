import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import math
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import distance_array as d_array
from MDAnalysis.analysis.contacts import contact_matrix
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF


struct_file = sys.argv[1]
traj_file = sys.argv[2]
cry_univ = mda.Universe(sys.argv[3],sys.argv[4])
time = sys.argv[5]

#Define MDAnalysis Universe
u = mda.Universe(struct_file,traj_file)
u_copy = u.copy()
#Protein selcetions
protein=u.select_atoms("byres name BB")     #All protein beads (1094)
protein_rms = u.select_atoms("byres name BB")

chainA_ref = cry_univ.select_atoms("bynum 52425-52882")
chainB_ref = cry_univ.select_atoms("bynum 52883-53340")
chainA = protein.select_atoms("bynum 52425-52882")
chainB = protein.select_atoms("bynum 52883-53340")
chainB_mic = protein_rms.select_atoms("bynum 52883-53340")

side_chains = protein.select_atoms("(byres name BB) and not name BB")  #Only side chain beads (590)
backbone = protein.select_atoms("name BB")    #Only backbone beads(504)
print(chainA)
#sys.exit()
#Select membrane atoms
memb = u.select_atoms("byres name PO4")     #All lipid beads (54728); Not a useful selection
#General parameters
po4 = u.select_atoms("name PO4")        #Selecting phospho head groups for all lipids (4414)
print(po4)
print(len(po4))
# sys.exit()
popc_lipid = u.select_atoms("resname POPC and name PO4") #Selecting phospho headgroups for diff lipids (3534)
pops_lipid = u.select_atoms("resname POPS and name PO4") #(440)
pop6_lipid = u.select_atoms("resname POP6 and name PO4") #(440)

com_po4 = po4.center_of_geometry()
upp_leaf_ind = []
low_leaf_ind = []
for i in range(len(po4)):
    z = po4[i].position[2]
    if z > com_po4[2]:
        upp_leaf_ind.append(i)
    elif z < com_po4[2]:
        low_leaf_ind.append(i)


upp_leaf = po4[upp_leaf_ind]
low_leaf = po4[low_leaf_ind]

# print(upp_leaf.positions[2454][2])
# print(low_leaf[0].position[2])
# sys.exit()
"""
Atom Group definitions for different patches
"""

"""
METAD CV patches
"""
#d1
#Chains A & B
#Sites 1 and 2 in each chain
#In this system, residue numbering of protein starts from 60797 (Atom 1 on chainA) and ends at  161872 (Atom 2 on chainB)
#Atom no. in this sytem is obtained by Atom no. in solution system - 9 + 60797  -> Chain A
#Atom no. in this sytem is obtained by Atom no. in solution system - 18 + 60797  -> Chain B

site_1A = protein.select_atoms("bynum 52436-52453") # First binding site in chain A         #Residues 57LYS-> 64THR  (pdb numbering)
site_1B = protein.select_atoms("bynum 52937-52968") #First binding site in chain B              #Residues 78GLU-92ASP
site_2A = protein.select_atoms("bynum 52799-52810") #Second binding site in chain A             #Residues 224LEU-229ASP
site_2B = protein.select_atoms("bynum 53178-53196") #Second binding site in chain B             #Residues 189ASN-196LYS

"""
Long axis along each chain
"""
#Chain A
#longA1 = protein.select_atoms("name BB and resid 180 and segid seg_0_Protein_A")
longA1 = protein.select_atoms("bynum 52668")
#longA1 = protein.select_atoms("bynum 245")
#longA2 = protein.select_atoms("name BB and resid 252 and segid A")
longA2 = protein.select_atoms("bynum 52819")
#longA2 = protein.select_atoms("name BB and resid 114 and segid seg_0_Protein_A")

#Chain B

#longB1 = protein.select_atoms("name BB and resid 432 and segid seg_1_Protein_B")
longB1 = protein.select_atoms("bynum 53126")
#longB1 = protein.select_atoms("bynum 714")
#longB2 = protein.select_atoms("name BB and resid 252 and segid B")
#longB2 = protein.select_atoms("name BB and resid 366 and segid seg_1_Protein_B")
longB2 = protein.select_atoms("bynum 53277")

"""
Membrane binding patches
"""
memb_A = protein.select_atoms("bynum 52433-52435 or bynum 52449-52451 or bynum 52456-52458 or bynum 52464-52466 or bynum 52594-52596 or bynum 52601-52603")
memb_B = protein.select_atoms("bynum 52891-52893 or bynum 52907-52909 or bynum 52914-52916 or bynum 52922-52924 or bynum 53052-53054 or bynum 53059-53061")

"""
Dimerisation Angle
"""
##Definning AtomGroup for Angle of Dimerisation
#theta1
#End atoms (two residues) of Helix 2 in each chain
hel_angle1A = protein.select_atoms("name BB and bynum 52430-52432")
hel_angle2A = protein.select_atoms("name BB and bynum 52657-52663")
hel_angle1B = protein.select_atoms("name BB and bynum 52888-52890")
hel_angle2B = protein.select_atoms("name BB and bynum 53115-53121")

#theta2
#Considering only the outer arms of helix2
out_hel_angle1A = protein.select_atoms("name BB and bynum 52613-52616")
out_hel_angle2A = protein.select_atoms("name BB and bynum 52657-52663")
out_hel_angle1B = protein.select_atoms("name BB and bynum 53071-53074")
out_hel_angle2B = protein.select_atoms("name BB and bynum 53115-53121")





"""
Restraining Wall patches
"""
# wall_A = protein.select_atoms("name BB and resid 5021-5024 and segid seg_620_Protein_P")
#wall_A = protein.select_atoms("name BB and resid 112-115")
# wall_B = protein.select_atoms("name BB and resid 5262-5265 and segid seg_621_Protein_Q")
#wall_B = protein.select_atoms("name BB and resid 364-367")


"""
Defining beads to calculate COM of each chain
"""
com_A_hi_flx = protein.select_atoms("bynum 52622 or bynum 52640-52643 or bynum 52660 or bynum 52680 or bynum 52696 or bynum 52833")
com_B_hi_flx = protein.select_atoms("bynum 52912 or bynum 53085 or bynum 53117 or bynum 53139 or bynum 53151 or bynum 53280 or bynum 53339")
#
com_A_low_flx = protein.select_atoms("bynum 52440 or bynum 52461 or bynum 52500 or bynum 52550 or bynum 52618 or bynum 52800 or bynum 52876")
com_B_low_flx = protein.select_atoms("bynum 52889 or bynum 52931 or bynum 53000 or bynum 53100 or bynum 53202 or bynum 53250 or bynum 53330")
#
com_A_med_flx = protein.select_atoms("bynum 52622 or bynum 52461 or bynum 52660 or bynum 52550 or bynum 52696 or bynum 52800")
com_B_med_flx = protein.select_atoms("bynum 52912 or bynum 52931 or bynum 53117 or bynum 53202 or bynum 53280 or bynum 53330")


"""
Defining beads for calculating pitch and yaw
"""
#Pitch
#Pick pair of beads to calculate the angle between the vector and z-axis
#Beads along the longitudnal axis



def check_image(pos1,pos2):
    del_x = pos1[0] - pos2[0]
    del_y = pos1[1] - pos2[1]
    del_z = pos1[2] - pos2[2]

    #print("Del_y:", del_y)
    #print("Del_z:", del_z)

    transVec = np.zeros((3))
    mic = False

    if abs(del_z) > box_dims[2]/2:
        mic = True
        #print("Z")
        if del_z > 0:
            transVec[2] = box_dims[2]
        else:
            transVec[2] = -box_dims[2]

    r_nopbc = ((del_x)**2 + (del_y)**2 + (del_z)**2)**0.5

    del_x = del_x - (box_dims[0])*round(del_x/box_dims[0],0)
    del_y = del_y - (box_dims[1])*round(del_y/box_dims[1],0)
    del_z = del_z - (box_dims[2])*round(del_z/box_dims[2],0)
    r = ((del_x)**2 + (del_y)**2 + (del_z)**2)**0.5

    return(r,r_nopbc,transVec,mic)

def calc_relVec(Ag1,com):
    pos_array = Ag1.positions
    return(pos_array - com)
# print(len(upp_leaf))
# print(len(low_leaf))

# print(len(memb))
# print(len(popc_lipid))
# print(len(pops_lipid))
# print(len(pop6_lipid))
# print(len(po4))
# print(len(protein))
# print(len(side_chains))
# print(len(backbone))
#sys.exit()
ts = u.trajectory.ts
global box_dims
box_dims = ts.dimensions
#Find COGs of each group
cog_protein = protein.center_of_geometry(compound='segments',unwrap=True)
cog_upp = upp_leaf.center_of_geometry()
cog_low = low_leaf.center_of_geometry()
cog_po4 = po4.center_of_geometry()

print(cog_protein)
print(cog_upp)
print(cog_low)
# sys.exit()
print(chainA.center_of_geometry())
print(chainB.center_of_geometry())


# u.transfer_to_memory()
n_frames = u.trajectory.n_frames
step = 200
chunks = np.arange(0,n_frames,step)
fr_count = 0
print(type(u.trajectory))
man_rmsf = []
chainA_rmsf = np.zeros(len(chainA))
chainB_rmsf = np.zeros(len(chainB))

# t1 = time_mod.perf_counter()
for chunk in chunks:
    u.transfer_to_memory(chunk,chunk+step)
    print("Analyzing frames: ",chunk , " - ", chunk+step)
    # t3 = time_mod.perf_counter()
    for i in range(chunk+step):
        protein.unwrap(reference='cog',inplace=True)

        comProA = chainA.center_of_geometry(pbc=False)
        comProB = chainB.center_of_geometry(pbc=False)

        new_comProA = comProA
        new_comProB = comProB


        # (d,d_nopbc,transVec,mic) = check_image(np.reshape(comProA,3),np.reshape(comProB,3))
        # relVec = calc_relVec(chainB,comProB)
        # new_comProB = comProB + transVec
        # new_positions_B = relVec + new_comProB
        # if mic:
        #     print("Before: ", protein_rms.positions[52882,:])
        #     itions_B
        #     chainB_mic.translate(transVec)
        #
        #     print("Positions after: ", protein_rms.positions[52882,:])
        #     print("After chainB_mic: ", chainB_mic.positions[0,:])

        # #Chain A RMSF
        # prealigner = align.alignto(u,u,select="name BB and bynum 52425-52882")
        # # reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
        # # reference = mda.Merge(protein).load_new(
        #             # reference_coordinates[:, None, :], order="afc")
        # aligner = align.alignto(u,cry_univ,select="name BB and bynum 52425-52882")
        #
        #
        # del_posA = np.sum((chainA.positions - chainA_ref.positions)**2)**0.5
        # chainA_rmsf = chainA_rmsf + del_posA**2
        # print("ChainB")
        # #Chain B RMSF
        # prealigner = align.alignto(u,u,select="name BB and bynum 52883-53340")
        # # reference_coordinates_B = u.trajectory.timeseries(asel=chainB).mean(axis=1)
        # # reference_B = mda.Merge(protein).load_new(
        # #             reference_coordinates_B[:, None, :], order="afc")
        # aligner = align.alignto(u,cry_univ,select="name BB and bynum 52883-53340")
        # del_posB = np.sum((chainB.positions - chainB_ref.positions)**2)**0.5
        # chainB_rmsf = chainB_rmsf + del_posB**2
        # print("Done")

chainA_rmsf = (chainA_rmsf/n_frames)**0.5
chainB_rmsf = (chainB_rmsf/n_frames)**0.5

protein.unwrap(reference='cog',inplace=True)

#Chain A RMSF
prealigner = align.AlignTraj(u,u,select="name BB and bynum 52425-52882",in_memory=True).run()
# reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
# reference = mda.Merge(protein).load_new(
            # reference_coordinates[:, None, :], order="afc")
aligner = align.AlignTraj(u,cry_univ,select="name BB and bynum 52425-52882", in_memory=True).run()
rmsf_A = RMSF(chainA,verbose=True).run()

#Chain B RMSF
prealigner = align.AlignTraj(u,u,select="name BB and bynum 52883-53340",in_memory=True).run()
# reference_coordinates_B = u.trajectory.timeseries(asel=chainB).mean(axis=1)
# reference_B = mda.Merge(protein).load_new(
#             reference_coordinates_B[:, None, :], order="afc")
aligner = align.AlignTraj(u,cry_univ,select="name BB and bynum 52883-53340", in_memory=True).run()
rmsf_B = RMSF(chainB,verbose=True).run()

# #TOtal chain RMSF
# prealigner = align.AlignTraj(u,u,select="name BB",in_memory=True).run()
# reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
# reference = mda.Merge(protein).load_new(
#             reference_coordinates[:, None, :], order="afc")
# aligner = align.AlignTraj(u,reference,select="name BB", in_memory=True).run()
# rmsf = RMSF(protein,verbose=True).run()

fig,ax = plt.subplots()
ax.plot(chainA.resnums,rmsf_A.rmsf,label='chainA')
ax.plot(chainB.resnums,rmsf_B.rmsf,label='chainB')
ax.plot(chainA.resnums,chainA_rmsf,label="ChainA-manual")
ax.plot(chainB.resnums,chainB_rmsf,label="ChainB-manual")
ax.legend()
plt.show()

#Identify regions of high flexibility
high_flexible_A = []
low_flexible_A = []
for i in range(len(rmsf_A.rmsf)):
    if rmsf_A.rmsf[i] <= 2.0:
        low_flexible_A.append(chainA[i].ix)
    elif rmsf_A.rmsf[i] > 4.5:
        high_flexible_A.append(chainA[i].ix)

high_flexible_B = []
low_flexible_B = []
for i in range(len(rmsf_B.rmsf)):
    if rmsf_B.rmsf[i] <= 2.0:
        low_flexible_B.append(chainB[i].ix)
    elif rmsf_B.rmsf[i] > 4.5:
        high_flexible_B.append(chainB[i].ix)

print("------- Chain A ------")
print("Atoms with high fluctuations: ", high_flexible_A)
print("Atoms with low fluctuations: ",low_flexible_A)

print("----- Chain B--------")

print("Atoms with high fluctuations: ", high_flexible_B)
print("Atoms with low fluctuations: ",low_flexible_B)
#The pickle file will store contacts data from each batch of trajectory
with open("RMSF_data_"+time,"w") as fl:
    fl.write("ResidueA\tRMSFA\tResidueB\tRMSFB\n")
    for i in range(len(rmsf_A.rmsf)):
        fl.write(str(chainA[i].ix))
        fl.write("\t")
        fl.write(str(rmsf_A.rmsf[i]))
        fl.write("\t")
        fl.write(str(chainB[i].ix))
        fl.write("\t")
        fl.write(str(rmsf_B.rmsf[i]))
        fl.write("\n")
