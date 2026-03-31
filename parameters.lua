-- ###########################################################################
-- 28-01-2026
-- Constants for ice sheet + bedrock setup (+ filler layer).
-- Parameter values chosen to correspond to Elisa's dimensionless values.
--
-- - Unit system = 
--     -- m
--     -- yr
--     -- MPa
--     -- K
-- ##########################################################################

-- # All units are coverted from SI to MPa, yr, m

-- # SI to MPa, yr, m conversions
yearindays = 365.25
yearinsec = yearindays*24*60*60
MPainPa = 1.0e6
GL=1.0526e+06

-- # Dimensionless parameters and prescribed scalings
G = 0.5   -- geothermal heat flux
alpha = 1  -- Brinkmann number (shear heating rate)
aND = 1      -- accumulation
Pe = 1     -- Peclet number (advection/diffusion)
TsND = -1    -- surface temperature
dbdx = 0.05    -- bed slope
Tscale = 10 --10  -- temperature scale (K)
gamma0 = 0.1   -- friction coefficient (0.1 fast, 3 slow)
delta = 0.03   -- subtemp sliding law exponent

-- # Physical parameters
rhoi = (920/MPainPa)/(yearinsec^2)    -- ice density (kg m^-3 = Pa m^2 s^2 --> MPa m^2 yr^2)
rhow = (1003.27/MPainPa)/(yearinsec^2)   -- water density (kg m^-3 = Pa m^2 s^2 --> MPa^ m^2 yr^2)
heatcapacity = 2097*yearinsec^2    -- specific heat capacity (J kg^-1 K^-1 = m s^-2 K^-1 --> m yr^-2 K^-1)
heatconductivity = 2.1/MPainPa*yearinsec   -- thermal conductivity (W m^-2 K^-1 = Pa m s^-1 K^-1 --> MPa m yr^-1 K^-1)
gravity = 9.81*yearinsec^2           -- acceleration due to gravity (m s^-2 --> m yr^-1)
viscosity = 1e14/MPainPa/yearinsec            -- viscosity of ice (Pa s --> MPa yr) - based on Nye 1969 value of 1 bar year = 3.15..e12 Pa s
T0 = 273.15               -- melting temperature (K)
Lw = 334000.0*yearinsec^2   -- latent heat of (melting) ice (J kg^-1 = m^2 s^-2 --> m^2 yr^-2)

-- # Calculated scalings and dimensional parameters
uscale = (alpha*heatconductivity*Tscale/viscosity)^(1/2)   -- velocity scale (m/yr)
xscale = rhoi*(heatcapacity^3)*(uscale^5)*(viscosity^2)/((Pe^3)*(heatconductivity^3)*(gravity^2))   -- x scale (m)
zscale = (viscosity*uscale*xscale/(rhoi*gravity))^(1/3)    -- z scale (m)

ascale = ((uscale^4)*viscosity/(rhoi*gravity*(xscale^2)))^(1/3);   -- accumulation scale (m yr^-1)
timescale = xscale/uscale
surfT = T0 - Tscale    -- surface temp (K)
gthf = G*heatconductivity*Tscale/zscale    -- geothermal heat flux (MPa m yr^-1)
aDIM = aND*ascale   -- accumulation (m yr^-1)    

-- # sliding law
dT = delta*Tscale    -- subtemp sliding exponent (K)
fcoeff = (gamma0*rhoi*gravity*zscale^2)/(uscale*xscale)    -- friction coefficient (MPa yr m^-1)

Source = 0.0

minh = 1.0   -- minimum ice thickness (m)

-- # parameters for GlaDS water sheet

p0 = 0.0   -- reference pressure (MPa)

-- # Bedrock
heatcapacityrock = 700*yearinsec^2    -- (J kg^-1 K^-1 = m s^-2 K^-1 --> m yr^-2 K^-1)) from toymodel.db in MLB example
heatconductivityrock = 3.0851/MPainPa*yearinsec    -- (J s^-1 m^-1 K^-1 = Pa m s^-1 K^-1 --> MPa m yr^-1 K^-1) from toymodel.db in MLB example
rhorock = (2780.0/MPainPa)/(yearinsec^2)   -- rock density (granite) (kg m^-3 = Pa m^2 s^2 --> MPa m^2 yr^2), taken from granite in the toymodel.db in the MLB example

Tsurf = T0 - Tscale

-- # Mesh Parameters
GL = 2.1*xscale
L = 2.5*xscale
GLRes =xscale/2000
BoundaryRes = xscale/30
nodes = 128



k = 0.01 -- sheet conductivity in m^(7/4) kg^(-1/2)
lr = 2 -- cavity spacing in m
hr = 0.1 -- Bedrock bump height in m

N0 = 1 -- Effective pressure constant (MPa) - GUESSED VALUE!
k0 = 4/3 -- Hydraulic conductivity exponent - value from Schoof and Mantelli 2021
D0 = 1e7 -- Hydraulic conductivity constant (??) - GUESSED VALUE!

Ncutepsdiff = 2/3 -- must be >1/3, used for regularisation
eps = 0.01 ---0.01 -- tolerance used in hydraulic conductivity parametrisation
Ncut = eps + Ncutepsdiff --0.001 -- tolerance used in hydraulic conductivity parametrisation
Ncoeff = (1/Ncut)*(1/3+((1/18)*(3*(Ncut-eps)-1))^(1/2))  -- for regularisation Ntilde = Ncoeff*N^2 + (1-2*delta*Ncoeff)*N + eps

