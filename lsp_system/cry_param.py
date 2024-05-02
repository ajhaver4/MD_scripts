import sys
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import math

struct_file = sys.argv[1]

#Define MDAnalysis Universe
u = mda.Universe(struct_file)

protein=u.select_atoms("byres name BB")


"""
Atom Group definitions for different patches
"""

"""
METAD CV patches
"""
#d1
#Chains A & B
#Sites 1 and 2 in each chain
site_1A = protein.select_atoms("bynum 12-29") # First binding site in chain A
site_1B = protein.select_atoms("bynum 524-556") #First binding site in chain B
site_2A = protein.select_atoms("bynum 383-394") #Second binding site in chain A
site_2B = protein.select_atoms("bynum 770-788") #Second binding site in chain B

"""
Long axis along each chain
"""
#Chain A
longA1 = protein.select_atoms("bynum 242")
#longA1 = protein.select_atoms("bynum 245")
#longA2 = protein.select_atoms("bynum 389")
longA2 = protein.select_atoms("bynum 92")

#Chain B

longB1 = protein.select_atoms("bynum 711")
#longB1 = protein.select_atoms("bynum 714")
#longB2 = protein.select_atoms("bynum 858")
longB2 = protein.select_atoms("bynum 561")

"""
Membrane binding patches
"""
memb_A = protein.select_atoms("bynum 9-11 or bynum 25-27 or bynum 32-34 or bynum 40-42 or bynum 165-167 or bynum 174-176 or bynum 181-183")
memb_B = protein.select_atoms("bynum 478-480 or bynum 494-496 or bynum 501-503 or bynum 509-511 or bynum 634-636 or bynum 643-644 or bynum 650-652")


def calc_dist(Ag1,Ag2,pbc):
    global box_dims
    s1A = Ag1.centroid(pbc=pbc,compound='segments')
    s1B = Ag2.centroid(pbc=pbc,compound='segments')


    del_x = s1A[0][0] - s1B[0][0]
    del_y = s1A[0][1] - s1B[0][1]
    del_z = s1A[0][2] - s1B[0][2]

    # del_x = del_x - (box_dims[0])*round(del_x/box_dims[0],0)
    # del_y = del_y - (box_dims[1])*round(del_y/box_dims[1],0)
    # del_z = del_z - (box_dims[2])*round(del_z/box_dims[2],0)


    r = ((del_x)**2 + (del_y)**2 + (del_z)**2)**0.5

    return(r)

def calc_angle(Ag1,Ag2,Ag3,Ag4):
    s1A = Ag1.centroid(pbc=False,compound='segments')
    s2A = Ag2.centroid(pbc=False,compound='segments')
    s1B = Ag3.centroid(pbc=False,compound='segments')
    s2B = Ag4.centroid(pbc=False,compound='segments')
    #print("Vector:",s1A)
    rA = s1A - s2A
    rB = s1B - s2B
    #print(rA)
    #print(rB)
    mod_rA = (rA[0][0]**2 + rA[0][1]**2 + rA[0][2]**2)**0.5
    mod_rB = (rB[0][0]**2 + rB[0][1]**2 + rB[0][2]**2)**0.5
    theta = math.acos((np.dot(rA.reshape(3),rB.reshape(3)))/(mod_rA*mod_rB))

    #dot_p = rA[0][0]*rB[0][0] + rA[0][1]*rB[0][1] +rA[0][2]*rB[0][2]
    #theta = math.acos((dot_p/(mod_rA*mod_rB)))

    return(theta)


d1 = calc_dist(site_1A,site_1B,False)
d2 = calc_dist(site_2A,site_2B,False)
l_angle = math.degrees(calc_angle(longA1,longA2,longB1,longB2))
dm = calc_dist(memb_A,memb_B,False)

print("----------------------------")
print("Crystal Structure Parameters")
print("MARTINI 3 model:")
print("----------------------------")
print("d1 = ",d1)
print("d2 = ",d2)
print("Axial angle = ", l_angle)
print("Dist between memb residues = ",dm)
