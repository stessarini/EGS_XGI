import numpy as np
import matplotlib.pyplot as plt
import string

def get_e_dep(MC_output_file_name):
    E_dep = np.array([])
    Std_Dev = np.array([])
    MC_output_file = open(MC_output_file_name, "r")
    line = MC_output_file.readline()
    found = False
    while(line and not found):
        line = line.strip()
        split_line = line.partition(" ")
        #print(split_line)
        if 'Deposited energy for' in line:
            #print(line)
            found = True
        line = MC_output_file.readline()
    line = MC_output_file.readline()
    is_edep_output = True
    while(line and is_edep_output):
        line = line.strip()
        split_line = line.partition(" ")
        #print(split_line)
        if split_line[0].isnumeric():
            #print(split_line[0])
            edep_line = split_line[2].strip()
            values = edep_line.partition('+/-')
            #print(values)
            e_dep = float(values[0].strip())
            std_dev = float(values[2].strip())
            #print(e_dep)
            #print(std_dev)
            E_dep = np.append(E_dep, e_dep)
            Std_Dev = np.append(Std_Dev, std_dev)
            #print(E_dep)
            #print(Std_Dev)
        else:
            is_edep_output = False
        line = MC_output_file.readline()
    return [E_dep, Std_Dev]

cases = string.ascii_uppercase[1:13:2]
print(cases)
E_dep = np.zeros((6, 9))
std_dev = np.zeros((6, 9))
#0: leave anywhere but through detector plane
#1: envelope
#2: wafer
#3: Air gaps (grating)
#4: periodic structure (grating) or all grating if there is no explicit grating structure
#5: Silicon cylinder
#6: Polystyrene cylinder
#7: exit through detector plane
#8: all the rest

#unfortunately have to assign the energies manually
name = '2mmx1mmFOV_PlaneWave_Case_'
#for i in cases[]:
#    file = name + i + '.log'
#    print(file)


#########################################
file = name + "B.log"
[Ed, sd] = get_e_dep(file)
E_dep[0,0] = Ed[0]
E_dep[0,1] = Ed[1]
E_dep[0,7] = Ed[2]
E_dep[0,8] = Ed[3]

std_dev[0,0] = sd[0]
std_dev[0,1] = sd[1]
std_dev[0,7] = sd[2]
std_dev[0,8] = sd[3]

#########################################
file = name + "D.log"
[Ed, sd] = get_e_dep(file)
E_dep[1,0] = Ed[0]
E_dep[1,1] = Ed[1]
E_dep[1,2] = Ed[2]
E_dep[1,7] = Ed[3]
E_dep[1,8] = Ed[4]

std_dev[1,0] = sd[0]
std_dev[1,1] = sd[1]
std_dev[1,2] = sd[2]
std_dev[1,7] = sd[3]
std_dev[1,8] = sd[4]

#########################################
file = name + "F.log"
[Ed, sd] = get_e_dep(file)
E_dep[2,0] = Ed[0]
E_dep[2,1] = Ed[1]
E_dep[2,2] = Ed[2]
#E_dep2,[7] = Ed[3]
E_dep[2,5] = Ed[4]
E_dep[2,6] = Ed[5]
E_dep[2,7] = Ed[6]
E_dep[2,8] = Ed[7]

std_dev[2,0] = sd[0]
std_dev[2,1] = sd[1]
std_dev[2,2] = sd[2]
#std_dev2,[7] = sd[3]
std_dev[2,5] = sd[4]
std_dev[2,6] = sd[5]
std_dev[2,7] = sd[6]
std_dev[2,8] = sd[7]

#########################################
file = name + "H.log"
[Ed, sd] = get_e_dep(file)
E_dep[3,0] = Ed[0]
E_dep[3,1] = Ed[1]
E_dep[3,2] = Ed[2]
E_dep[3,4] = Ed[7]
E_dep[3,5] = Ed[4]
E_dep[3,6] = Ed[5]
E_dep[3,7] = Ed[6]
#E_dep[8] = Ed[7]

std_dev[3,0] = sd[0]
std_dev[3,1] = sd[1]
std_dev[3,2] = sd[2]
std_dev[3,4] = sd[7]
std_dev[3,5] = sd[4]
std_dev[3,6] = sd[5]
std_dev[3,7] = sd[6]
#std_dev[8] = sd[7]


#########################################
file = name + "J.log"
[Ed, sd] = get_e_dep(file)
E_dep[4,0] = Ed[0]
E_dep[4,1] = Ed[1]
E_dep[4,2] = np.sum(Ed[2:2012])
E_dep[4,3] = np.sum(Ed[2012:4021:2])
E_dep[4,4] = np.sum(Ed[2013:4022:2])
E_dep[4,5] = Ed[4022]
E_dep[4,6] = Ed[4023]
E_dep[4,7] = Ed[-2]
E_dep[4,8] = Ed[-1]

std_dev[4,0] = sd[0]
std_dev[4,1] = sd[1]
std_dev[4,2] = np.sum(sd[2:2012])
std_dev[4,3] = np.sum(sd[2012:4021:2])
std_dev[4,4] = np.sum(sd[2013:4022:2])
std_dev[4,5] = sd[4022]
std_dev[4,6] = sd[4023]
std_dev[4,7] = sd[-2]
std_dev[4,8] = sd[-1]


#########################################
file = name + "L.log"
[Ed, sd] = get_e_dep(file)
E_dep[5,0] = Ed[0]
E_dep[5,1] = Ed[1]
E_dep[5,2] = np.sum(Ed[2:2012])
E_dep[5,3] = np.sum(Ed[2012:4021:2])
E_dep[5,4] = np.sum(Ed[2013:4022:2])
E_dep[5,5] = Ed[4022]
E_dep[5,6] = Ed[4023]
E_dep[5,7] = Ed[-2]
E_dep[5,8] = Ed[-1]

std_dev[5,0] = sd[0]
std_dev[5,1] = sd[1]
std_dev[5,2] = np.sum(sd[2:2012])
std_dev[5,3] = np.sum(sd[2012:4021:2])
std_dev[5,4] = np.sum(sd[2013:4022:2])
std_dev[5,5] = sd[4022]
std_dev[5,6] = sd[4023]
std_dev[5,7] = sd[-2]
std_dev[5,8] = sd[-1]

print(E_dep)
print(std_dev)
