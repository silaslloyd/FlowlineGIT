import numpy as np
from collections import defaultdict
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import glob
import os

def read_partitioned_mesh_nodes(base_dir):
    nodes = {}

    node_files = sorted(glob.glob(os.path.join(base_dir, "part.*.nodes")))
    if not node_files:
        raise RuntimeError("No part.*.nodes files found")
    print(node_files)

    for fname in node_files:
        L = input()
        with open(fname) as f:
            for line in f:
                print(line)
                if not line.strip():
                    continue
                parts = line.split()    #Node line format: Node ID : Partition index : X : Y : Z
                nid = int(parts[0])
                print(nid)
                x, y, z = map(float, parts[2:5])
                nodes[nid] = np.array([x, y, z])
    plt.pause(1000)
    return nodes

def collect_global_boundary_nodes(partition_dir):
    boundary_nodes = defaultdict(set)

    for fname in sorted(glob.glob(os.path.join(partition_dir, "part.*.boundary"))):
        with open(fname) as f:
            for line in f:
                parts = line.split()    #Boundary Line format:  Boundary Element ID : Boundary Number : Parent Element 1 : Parent Element 2 : Element Type : Element Nodes
                bid = int(parts[1])
                for nid in parts[5:]:
                    boundary_nodes[bid].add(int(nid))   #Nodes corresponding to boundary element

    if not boundary_nodes:
        raise RuntimeError("No boundary data found")

    return boundary_nodes


def find_free_surface_boundary(boundary_nodes, zs_dof):
    bid = min(
        boundary_nodes.keys(),
        key=lambda k: abs(len(boundary_nodes[k]) - zs_dof)
    )
    print(bid)
    return bid, boundary_nodes[bid]

def read_partitioned_zs(result_base, zs_name="zs"):
    zs_global = []
    permtable = []

    result_files = sorted(glob.glob(result_base + ".*"))
    if not result_files:
        raise RuntimeError("No partitioned result files found")

    
    for fname in result_files:
        reading = False
        with open(fname) as f:
            for line in f:
                if line.strip() == zs_name:  
                    reading = True
                    continue
                if reading:
                    parts = line.split()
                    if line.strip() == "zs residual":
                        break  
                    if len(parts) == 2:
                        permtable.append([parts[0],parts[1]])
                    if len(parts) == 1:
                        zs_global.append(float(parts[0]))

    permtable = np.array(permtable)
    idx =np.int16(permtable[:,1])
    zs_global = np.array(zs_global)
    zs_global = zs_global[idx]
    permtable[:,1] = zs_global

    return np.array(permtable)


def read_partitioned_zb(result_base, zb_name="zb"):
    zb_global = []
    permtable = []

    result_files = sorted(glob.glob(result_base + ".*"))
    if not result_files:
        raise RuntimeError("No partitioned result files found")

    for fname in result_files:
        reading = False
        with open(fname) as f:
            for line in f:
                if line.strip() == zb_name:  
                    reading = True
                    continue
                if reading:
                    parts = line.split()
                    if line.strip() == "zb residual":
                        break  
                    if len(parts) == 2:
                        permtable.append([parts[0],parts[1]])
                    if len(parts) == 1:
                        zb_global.append(float(parts[0]))

    permtable = np.array(permtable)
    idx =np.int16(permtable[:,1])
    zb_global = np.array(zb_global)
    zb_global = zb_global[idx]
    permtable[:,1] = zb_global

    return np.array(permtable)

def build_surface_interpolator(nodes, surface_nodes, zs):
    # Coordinates of surface nodes
    print(nodes)
    print(zs)
    plt.pause(1000)
    vals = []
    zvals = []
    surface_nodes = zs[:,0]

    for nid in surface_nodes:
        try:
            vals.append(nodes[nid][0])
        except KeyError:
            pass
    x = np.array(vals)

    for nid in surface_nodes:

        for row in zs:
            if row[0] == nid:
                zvals.append(row[1])     
    z= np.array(zvals)

    # Sort surface nodes spatially
    idx = np.argsort(x)

    x_sorted = x[idx]
    z_sorted = z[idx]
    print(np.shape(x_sorted))
    print(np.shape(z_sorted))


    z_of_x = interp1d(
        x_sorted,
        z_sorted,
        kind="linear",
        fill_value="extrapolate"
    )

    return x_sorted, z_sorted, z_of_x

# -------- PATHS --------
PARTITION_DIR = "/import/ontap-m-glaciology/Lloyd/Flowline/TillyElmer/footprint/partitioning.16"
RESULT_BASE   ="/import/ontap-m-glaciology/Lloyd/Flowline/TillyElmer/footprint/geometry_Stokes_heat.result"


ZS_DOF = 5998

# Load mesh nodes
# Load global mesh
nodes = read_partitioned_mesh_nodes(PARTITION_DIR)

# Collect boundaries
boundary_nodes = collect_global_boundary_nodes(PARTITION_DIR)
print(boundary_nodes)
# Identify free surface
surface_id, surface_nodes = find_free_surface_boundary(
    boundary_nodes, ZS_DOF
)
print(len(surface_nodes))
print(f"Free surface boundary ID: {surface_id}")
print(f"Surface nodes (unique): {len(surface_nodes)}")

# Read zs
# Read partitioned zs
zs = read_partitioned_zs(RESULT_BASE)
zb= read_partitioned_zb(RESULT_BASE)

print(f"Total zs DOFs read: {len(zs)}")


# Build interpolator
# Build interpolator
x_surf, z_surf, z_of_x = build_surface_interpolator(
    nodes, list(surface_nodes), zs
)
x_surf, z_base, z_of_x = build_surface_interpolator(
    nodes, list(surface_nodes), zb
)
plt.figure(figsize=(8, 4))
plt.scatter(x_surf, z_surf, s=6, label="Elmer free surface")
plt.scatter(x_surf, z_base, s=6, label="Elmer free surface")

xq = np.linspace(x_surf.min(), x_surf.max(), 1000)
#lt.plot(xq, z_of_x(xq), "r", lw=2, label="Interpolated")
plt.xlabel("x")
plt.ylabel("z_s")
plt.legend()
plt.tight_layout()
plt.show()
