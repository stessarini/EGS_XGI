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
