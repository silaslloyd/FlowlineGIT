import numpy as np
from collections import defaultdict
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import glob
import os


def read_partitioned_vars(result_base):
    data = defaultdict(lambda: defaultdict(list))
    VarNames = []
    Timesteps = []
    permtable = []
    current_vals = []
    currentvar = []
    result_files = sorted(glob.glob(result_base + ".*"))
    if not result_files:
        raise RuntimeError("No partitioned result files found")

    for fname in result_files:
        
        #if fname != "/import/ontap-m-glaciology/Lloyd/Flowline/FlowlineScaled/footprint/InitialT2Full.result.10" and fname != "/import/ontap-m-glaciology/Lloyd/Flowline/FlowlineScaled/footprint/InitialT2Full.result.11":
        #    continue
        print(fname)
        get_var_list = False
        get_var_vals = False
        read_perm_table = False
    
        with open(fname) as f:
            for line in f:

                if line.strip() == "Degrees of freedom:":
                    get_var_list = True
                    continue
                if line.split()[0] == "Total":
                    get_var_list = False
                if get_var_list == True:
                    parts = line.split()
                    VarString = ' '.join(parts[:parts.index(':')])
                    if VarString not in VarNames:
                        VarNames.append(VarString)

                if line.split()[0] == "Time:":
                    currentTime = line.split()[3]
                    if currentTime not in Timesteps:
                        Timesteps.append(currentTime)

                    get_var_vals = True
                    continue
                if get_var_vals == True:
     
                    if line.strip() in VarNames and line.strip() != currentvar and currentvar in VarNames:
                        idx = np.int16(np.array(permtable)[:,1])-1

                        current_vals = np.array(current_vals)[idx]

                       # if currentvar == "nodalcoords 1":
                       #     print(np.shape(current_vals))
                        nodeids = np.array(permtable)[:,0]
                        current_vals = np.column_stack((nodeids, current_vals))

                       # if currentvar == "nodalcoords 1":
                       #     print("CurrentValsNodes")
                       #     print(np.shape(current_vals))

                       #     print(np.shape(data[currentvar][currentTime]))
                        #print(data[currentvar])
                        if currentTime not in data[currentvar]:
                            data[currentvar][currentTime] = current_vals
                        else:
                            data[currentvar][currentTime] = np.vstack((data[currentvar][currentTime], current_vals))

                       # if currentvar == "nodalcoords 1":
                       #     print(np.shape(data[currentvar][currentTime]))

                        current_vals = []
                        currentvar = line.strip()
                    
                    if line.strip() in VarNames and currentvar not in VarNames:
                        currentvar = line.strip()
                        continue

                    if line.split()[0] == "Perm:" and line.split()[1] != "use":     
                        permtable = []   
                        read_perm_table = True
                        continue
                    if len(line.split()) == 1 and line.strip() not in VarNames:
                        read_perm_table = False
                        current_vals.append(np.float64(line.split()[0]))

                    if read_perm_table == True:
                        permtable.append([line.split()[0], line.split()[1]])
    print(VarNames) 
    return data, Timesteps


# -------- PATHS --------
RESULT_BASE   ="/import/ontap-m-glaciology/Lloyd/Flowline/FlowlineScaled/footprint/InitialT2Full.result"

data, Timesteps = read_partitioned_vars(RESULT_BASE)

for var in data:
    print(var)
    for t in data[var]:
        print("   ", t)

for T in Timesteps:

    zs = np.float64(data["zs"][T])
    xs = np.float64(data["nodalcoords 1"][T])
    ys = np.float64(data["nodalcoords 2"][T])
    zs = np.float64(data["nodalcoords 3"][T])
    Temperature = np.float64(data["temperature"][T])

    x_lookup = {nid: x for nid, x in xs}
    z_lookup = {nid: z for nid, z in zs}

    common_nodes = sorted(set(x_lookup) & set(z_lookup))

    x = np.array([x_lookup[n] for n in common_nodes])
    z = np.array([z_lookup[n] for n in common_nodes])
    print(Temperature)
    T_lookup = {nid: v for nid, v in Temperature}
    print(common_nodes)
    print(T_lookup)
    f = np.array([T_lookup[n] for n in common_nodes])

    plt.tricontourf(x, z, f, levels=20)
    plt.colorbar()
    plt.xlabel("x")
    plt.ylabel("z")
    plt.show()

plt.show()