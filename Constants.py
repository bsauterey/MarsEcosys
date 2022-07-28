"""
This code is under MIT license. See the License.txt file.
Environmental constants

Boris Sauterey
boris.sauterey@ens.fr
"""

import numpy as np 

###Environmental constants
T          = 253                   # Atmospheric Temperature (K)
Av         = 6e23                  # Avogadro number
p          = 1e5                   # Atmospheric pressure (Pa)     
R          = 8.31                  # Perfect gaz constant J.(K.mol)-1
Toc        = T                     # Oceanic Temperature (K)
S          = 1.4e18                # Ocean surface (cm2)
Voc        = 3e21                  # Ocean Volume (L)
g          = 3.711                 # Gravitational constant (m.s-2)
Mair       = 44e-3                 # Air molar mass (kg.mol-1)
M_CO2      = 44e-3                 # CO2 molar mass (kg.mol-1)
M_N2       = 28e-3                 # N2 molar mass (kg.mol-1)
M_CH4      = 16e-3                 # CH4 molar mass (kg.mol-1)
M_H2       = 2e-3                  # H2 molar mass (kg.mol-1)
rho        = p/((R/Mair)*T)        # Air density (kg.m^-3)
Cair       = rho/Mair*1e-3         # Mol of air per L-1
#Cair       = rho/Mair*1e-6        # Mol of air per cm-3
ntot       = (p*S*1e-4)/(g*Mair)   # Total number of air mol
a          = 1.438e15              # Intercept of the Photolysis of CH4 (yr-1)
b          = 2.0291                # Exponent of the Photolysis of CH4
Vco        = 0.63                  # Speed of deposition of CO (cm.yr-1) 
E          = 13.1e-4               # Escape rate of H2 to space (mol.cm-2.yr-1)
Volc       = 3*10**(10)/Av*365*24*60*60  # Volcanism (mol.cm-2.yr-1); pH2 = 800 ppm
Rec_prop   = 0.98                  # Recycling proportion
tau_oc     = 5000                  # Convective constant (yr)
totZ       = Voc*1e3 / S / 1e2     # Total ocean depth (m)
MZZ        = 100                   # Mixing zone depth
DOZ        = totZ - MZZ            # Deep ocean depth
Fconv      = (MZZ+DOZ)/tau_oc      # Mixing constant MZ <-> DO (m.yr-1)  

### Physical properties of elements
# Molar mass (g.mol-1)
MH         = 1
MC         = 12
MG         = 16
MCO        = 28
MN2        = 28
MCH3COOH   = 60

# Piston velocities (cm.s-1 -- assuming a 40µm mixing layer)
vH         = 1.3e-2                # H2
vC         = 4.8e-3                # CO2
vG         = 4.5e-3                # CH4
vCO        = 4.8e-3                # CO
vN2        = 4.8e-3                # CO

# Raining out of NOxs (mol.cm-2.s-1 (Wong et al 2017))
#vNO3       = 7.4e7/Av
#vNO2       = 1.4e7/Av
#vNO3       = 2.5e6/Av
#vNO2       = 3e4/Av
#vNO3       = 1.5e5/Av
#vNO2       = 1.5e3/Av
vNO3       = 9.3e4/Av
vNO2       = 7.8e1/Av

# Raining out of H2SO4 (mol.cm-2.s-1 (Ono et al 2003))
vH2SO4     = 1e9/Av

# Flux of water through hydrothermal vents (eliminates NOx through serpentinization d-1 (Wong et al 2017)) 
Fhyd       = 7.2e12/Voc/365 

# Fluxes (d-1 -- assuming mixing in 1x1x10000 cm water column)
QH         = vH   * 60*60*24 * 1e-4
QC         = vC   * 60*60*24 * 1e-4
QG         = vG   * 60*60*24 * 1e-4
QCO        = vCO  * 60*60*24 * 1e-4
QN         = 1e-5                  
QCH3COOH   = 0

# Fluxes (mol.L-1.d-1 -- assuming mixing in 1x1x10000 cm water column)
QNO3       = vNO3   * 60*60*24 * 1e-1
QNO2       = vNO2   * 60*60*24 * 1e-1
QH2SO4     = vH2SO4 * 60*60*24 * 1e-1
QN2        = vN2    * 60*60*24 * 1e-1

# Solubilities (mol.L-1.bar-1)
alphaH     = 7.8e-4
alphaG     = 1.4e-3
alphaC     = np.exp(9345.17/T-167.8108+23.3585*np.log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35.0)
alphaCO    = 1e-3
alphaN2    = 7e-4

### Prebiotic atmospheric composition (ppm)
#Gatm       = ntot*1e-20 *1e-6
Gatm       = ntot*1e3   *1e-6
#Catm       = ntot*1e3  *1e-6
Catm       = ntot*2500  *1e-6
pCO2       = Catm/ntot  *1e6
COatm      = ntot*1e1   *1e-6
N2atm      = ntot*9e5   *1e-6

### Prebiotic NOx composition (mol.L-1 (Wong et al 2017))
NO3inf     = vNO3*60*60*24*S/Voc /(Fhyd) 
NO2inf     = vNO2*60*60*24*S/Voc /(Fhyd)
H2SO4inf   = vH2SO4*60*60*24*S/Voc /(Fhyd) 


######################### Section that covers the C cycle ##################

R = 8.314 
T_0 = T	#
J = 7e15 	# diffusion coefficient between ocean and pore in kg/yr
zeta=0.3 	#CO2-dependence weathering
nout=0.35	#exponent for outgassing
mout=1.5	#exponent for outgassing
beta=0.1	#exponent for spreading rate
gam_ma=0.25	#H+-dependence oceanic basalt weathering
Te= 25.0 		#temperatude-dependence silicate weathering in K
K=77.8  	#conductivity (m/K)
Ebas=80.0e3		#energy activation J/mol
n=1.7		#exponent for carbonate precipitation

Mo=Voc 	#mass of the ocean
Mp=1.35e19 	#mass of pore
so=ntot/1e5/Mo 	# moles/pa/kg
tini=3.8	#Ga
#tini=2.6	#Ga
#tini=2.3	#Ga
#tini=3	#Ga
tgrowth =2.5 	#Time of continental growth in Ga
tbio=0.6	#Time of biological weathering in Ga
Bprecambrian=0.5 #Weathering factor at 4 Ga
fsed=0.6	#relative depth of sediment at 4 Ga
Larchean=0.2	#fraction emerged land during Archean compared to today
concHo0_hadean=10**(-6.5)
Fout0 = 8.0e12 	#mol/yr present volcanism
Fcarb0 =10.5e12 #mol/yr present carbonate weathering
Fsil0 =10.5e12  #mol/yr present silicate weathering
Fdiss0 = 0.45e12 #mol/yr present seafloor weathering
Tpore0 = 367.0	#present temperature of pore in K
Pocean0= 0.45e12 #mol/yr present ocean precipitation
Ppore0=	0.45e12	#mol/yr present pore precipitation
#omegap0=	#
#omegao0=	#
#concCaini=10e-3 #concentration Ca2+ (mol/L)
concHo0=6e-9	#present [H+] mol/L
#omegap0=0	#
#omegao0=0	#