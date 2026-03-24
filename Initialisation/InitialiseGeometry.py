import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import re
import math
from scipy.optimize import root_scalar


def load_lua_constants(filepath):
    variables = {}
    
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        # Remove comments
        line = line.split("--")[0].strip()
        if not line:
            continue
        
        # Skip print statements
        if line.startswith("print"):
            continue
        
        # Replace Lua power operator ^ with Python **
        line = line.replace("^", "**")
        
        # Replace Lua math syntax if needed
        line = line.replace(";", "")
        
        # Execute in controlled namespace
        try:
            exec(line, {"__builtins__": None, "math": math}, variables)
        except Exception:
            pass  # skip anything unexpected
    
    return variables


Plot = False
### Load Scaling Parameters from parameters.lua ###

params = load_lua_constants("../parameters.lua")
viscosity = params["viscosity"]
rhoi = params["rhoi"]
gravity = params["gravity"]
rhow = params["rhow"]
C = params["fcoeff"]
a = params["aDIM"]
Tscale = params["Tscale"]
T0 = params["T0"]


# Load Dimensionless SIA steady state geometry
x = np.loadtxt("Dimensionless_ICs/x_SIA.txt")
h = np.loadtxt("Dimensionless_ICs/h_SIA.txt")
T = np.loadtxt("Dimensionless_ICs/T_SIA.txt")
h = (np.nan_to_num(h, nan=0))
T = (np.nan_to_num(T, nan=0))
h = np.maximum(h, 1e-5)

bedx = np.loadtxt("Dimensionless_ICs/bed_x.txt")
bedz = np.loadtxt("Dimensionless_ICs/bed_z.txt")

# Create functional form for scaled bedrock
bed = lambda x:(params["zscale"]*(bedz[1] - bedz[0]))/(params["xscale"]*(bedx[1] - bedx[0]))*x

# Scale sheet geometry
h = params["zscale"]*h
x = params["xscale"]*x
T = np.minimum(T*Tscale + T0,T0)

#Functional form for sheet geometry
h_func = CubicSpline(x, h)

#Functional form for temperature
T_func = CubicSpline(x, T)

#Domain Size for function
x_smaple = np.linspace(0,2.1,100000)*params["xscale"]

#Calcualte steady state grounding line through Schoof condition
coeff = (((rhoi*gravity)**2 * (1-rhoi/rhow))/(8*C*viscosity))**(0.5) 
f = lambda x: coeff*h_func(x)**(5/2) - a*x
res = root_scalar(f, bracket=[0, 3*params["xscale"]], method='brentq')
xg = res.root

#Calculate required bedrock height at grounding line through flotation condition and shift bedrock function
b_xg = - rhoi/rhow * h_func(xg)
bed_shifted = lambda x: bed(x) - bed(xg) + b_xg

#Calculate analytic extensional shelf
gamma = rhoi*(1-rhoi/rhow)*gravity/(8*viscosity*a)
xshelf = np.linspace(xg, 5*params["xscale"], 100000)
F = (xshelf/xg)**2 * h_func(xg)**2 / (1-gamma*h_func(xg)**2)
hshelf = np.sqrt(np.abs(F/(1+gamma*F)))

#Create functional forms for shelf upper and lower shelf surface
s_shelf = CubicSpline(xshelf,hshelf*(1-rhoi/rhow)*(xshelf>xg))
l_shelf  = CubicSpline(xshelf,hshelf*(-rhoi/rhow))

#Create functional forms for upper and lower grounded surface
s_grounded = lambda x: h_func(x) + bed_shifted(x)
l_grounded = lambda x: bed_shifted(x)

#Create full functional forms for upper and lower ice surfaces
s = lambda x: np.maximum(s_shelf(x), s_grounded(x))
l = lambda x: np.maximum(l_shelf(x), l_grounded(x))


# Plot

X = np.linspace(0,6E6,100000)
zs = np.column_stack((X, s(X)))
zb = np.column_stack((X, l(X)))
b = np.column_stack((x, bed_shifted(x)))
T = np.column_stack((X, T_func(X)))

np.savetxt('Dimensional_ICs/zs_dim.xy', zs, fmt='%g', delimiter=' ')
np.savetxt('Dimensional_ICs/zb_dim.xy', zb, fmt='%g', delimiter=' ')
np.savetxt('Dimensional_ICs/bed_dim.xy', b, fmt='%g', delimiter=' ')
np.savetxt('Dimensional_ICs/Tb_dim.xy', T, fmt='%g', delimiter=' ')

if Plot == True:
    plt.figure()
    var1 = np.loadtxt("Dimensional_ICs/zs_dim.xy")
    var2 = np.loadtxt("Dimensional_ICs/Tb_dim.xy")

    plt.plot(var1[:,0], var1[:,1])
    plt.plot(var2[:,0], var2[:,1])
    plt.show()