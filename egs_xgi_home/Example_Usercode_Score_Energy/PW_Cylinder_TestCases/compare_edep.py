import numpy as np
import matplotlib.pyplot as plt
from get_e_dep import get_e_dep



name = '2mmx1mmFOV_PlaneWave_Case_'
E_dep = np.zeros((6, 9))
std_dev = np.zeros((6, 9))

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


np.save("E_dep_GI.npy",E_dep)
np.save("std_dev_GI.npy",std_dev)
