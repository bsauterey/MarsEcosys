from Bio import *
from Constants import *
from Climate import *
from Photochem import *
from math import *
import gzip
import pickle
import time as libtime
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess as sp
from stepfunc import *
from scipy.integrate import solve_ivp
import sys
from T_to_rho_and_Tsurf import *


## Sets environmental variables and parameters
# Atmosphere
p = np.random.uniform(0.5,3)*10**5
pCH4           = 1e-4
pH2            = np.random.uniform(3e-3,0.1)-2*pCH4
pCO            = 1e-10
pCO2           = (1-(pCH4+pH2+pCO))*0.95
pN2            = (1-(pCH4+pH2+pCO))*0.05
ntot           = (p*S*1e-4)/(g*(pCO2*M_CO2+pN2*M_N2+pCH4*M_CH4+pH2*M_H2))
Catm           = pCO2 * ntot * 1e-6
Hatm           = pH2  * ntot * 1e-6
Gatm           = pCH4 * ntot * 1e-6
N2atm          = pN2  * ntot * 1e-6
COatm          = pCO  * ntot * 1e-6

# Brines' freezing point
T_freez        = 252    #Can be 203, 252, or 273

# Initial average T°, ice-free portion of the surface, and T° in that region
T_av           = T_func(p*1e-5,pH2,pCH4)
rho            = float(T_to_rho252(T_av)) #Can be T_to_rho203, T_to_rho252, or T_to_rho273
T_surf         = max(T_av,float(T_to_Tsurf252(T_av))) #Can be T_to_Tsurf203, T_to_Tsurf252, or T_to_Tsurf273

# Soil properties
a_eps          = np.random.uniform(0.2,0.6)
a_to           = np.random.uniform(1.5,2.5)
r_surf         = np.exp(np.random.uniform(log(1e-4),log(1e-3)))
ztot           = np.random.uniform(5,20)
layers         = 100
a_T            = np.random.uniform(10,30)
z              = ztot
Ti             = T_surf + a_T * ztot
z_vec          = np.linspace(0,ztot,layers)
T_vec          = T_surf + a_T * z_vec
par            = [a_eps,a_to,r_surf]

# Fluxes in and out of the atmosphere
E              = 1.6e13
Esc            = E*(pH2+pCH4/2)
F_CH4          = float(FCH4_func(pH2,pCH4))
F_H2           = float(FH2_func(pH2,pCH4))
Volc           = Esc - F_H2

# Vector containing the environmental variables
Env            = [pH2*p*1e-5*alphaH*1e-3*Av,pCO2*p*1e-5*np.exp(9345.17/T_surf-167.8108+23.3585*np.log(T_surf)+(0.023517-2.3656e-4*T_surf+4.7036e-7*T_surf**2)*35.0)*1e-3*Av,pCH4*p*1e-5*alphaG*1e-3*Av,pN2*p*1e-5*alphaN2*1e-3*Av]

## Sets the initial conditions and shape them for the numerical solver
# Environmental variables
init_chem      = []
for i in np.arange(len(Env)):
	init_chem += [Env[i] for j in np.arange(layers)]

# Biological variables
FTnames        = ['meth','NO3metht','NO2metht','H2SO4metht','acet','acet2','acett','ferm','photoH2']
FTlist         = {'meth':True,'NO3metht':False,'NO2metht':False,'H2SO4metht':False,'acet':False,'acet2':False,'acett':False,'ferm':False,'photoH2':False}
t_intro_list   = {'meth':1e100,'NO3metht':1e100,'NO2metht':1e100,'H2SO4metht':1e100,'acet':1e100,'acet2':1e100,'acett':1e100,'ferm':1e100,'photoH2':1e100}
traits_list   = {'meth':[],'NO3metht':[],'NO2metht':[],'H2SO4metht':[],'acet':[],'acet2': [],'acett':[],'ferm':[],
                 'photoH2':[]}

for i in FTnames:
	if FTlist[i] is True :
		traits_list[i] = Run_def_FT(i,T_vec)   ### SHOULD BE THE DEFAULT !!!
		#traits_list[i] = Run_def_FT(i,T_vec,Forced_rc=True)

init_NC    = []
init_X0    = []
for i in FTnames:
	if FTlist[i] is True:
		init_NC = [1e-99 for j in np.arange(layers)]
		for j in np.arange(layers):
			traits   = Run_def_FT(i,T_vec[j])   ### SHOULD BE THE DEFAULT !!!
			#traits   = Run_def_FT(i,T_vec[j],Forced_rc=True)
			init_X0 += [traits[9]*10]

# Total vector (environment then cell count then cell biomass content)
init = np.array(init_chem + init_NC + init_X0)

## Saves the values of the main parameters of the simulation under the name that must be provided as argument
pickle.dump([p,pH2,pCO2,pCH4,a_eps,a_to,r_surf,T_av,ztot,a_T],open(sys.argv[1]+'.params','wb'))

der = stepfunc(0,0,init,layers,FTnames,FTlist,traits_list,Env,z_vec,T_vec,par,ztot,p)

X0T     = [[] for i in np.arange(layers)]
NCT     = [[] for i in np.arange(layers)]
HT      = [[] for i in np.arange(layers)]
CT      = [[] for i in np.arange(layers)]
NT      = [[] for i in np.arange(layers)]
GT      = [[] for i in np.arange(layers)]

y       = [init]
y1      = init
Env1    = Env
ntot1   = ntot
p1      = p
T_vec1  = T_vec
T_surf1 = T_surf
T_av1   = T_av
rho1    = rho
pH2_1   = pH2
pCH4_1  = pCH4
pCO2_1  = pCO2
pN2_1   = pN2
CH4tot1 = ntot*pCH4_1
H2tot1  = ntot*pH2_1

pT      = [p1]
T_surfT = [T_surf1]
T_avT   = [T_av1]
rhoT    = [rho1]
pCH4T   = [pCH4_1]
pH2T    = [pH2_1]
pCO2T   = [pCO2_1]
pN2T    = [pN2_1]
DiffHT  = [0]
DiffGT  = [0]
t       = [0]

t_count = 0
ratchet_step = 1000          # Determines how frequently the state of the system is saved
ratchet = ratchet_step

#tmax = 365*300
tmax = 365*400
#tmax = 365*1
#tmax = 10000
if rho1 > 0:                 # Check whether Mars is initially habitable or not
	while t[-1] < tmax:      # Time loop

		y0        = y1
		Env0       = Env1
		p0        = p1
		ntot0     = ntot1
		T_vec0    = T_vec1
		T_surf0   = T_surf1
		T_av0     = T_av1
		rho0      = rho1
		pCH4_mol0 = Env0[2]
		pH2_mol0  = Env0[0]
		pH2_0     = pH2_1
		pCH4_0    = pCH4_1
		pCO2_0    = pCO2_1
		pN2_0     = pN2_1
		CH4tot0   = CH4tot1
		H2tot0    = H2tot1

		for i in FTnames:
			if FTlist[i] is True :
				traits_list[i] = Run_def_FT(i,T_vec0)

		# Evaluates soil ecosystem derivatives, then new state
		## Derivative
		der = stepfunc(t_count,ratchet,y0,layers,FTnames,FTlist,traits_list,Env0,z_vec,T_vec0,par,ztot,p0)
		no_zero  = np.where(np.array(der) != 0)
		if t[-1] < 2000:
			dt   = max(min(abs(y0[no_zero]/np.array(der)[no_zero]))/10,1e-2)    # Adaptive time step (must be low at first)
		else:
			#dt   = min(abs(y0[no_zero]/np.array(der)[no_zero]))/10
			dt   = max(min(abs(y0[no_zero]/np.array(der)[no_zero]))/5,1e-1)     # Adaptive time step
		if t[-1] + dt > tmax: dt = tmax - t[-1]
		## New state
		y1             = y0 + np.array(der) * dt
		pCH4_mol_surf1 = y1[(2*layers):(2*layers+layers)][0]
		pH2_mol_surf1  = y1[(0*layers):(0*layers+layers)][0]

		# Evaluate the atmosphere derivatives, then new state
		## Photochemistry
		F_CH4 = float(FCH4_func(pH2_0,pCH4_0))
		F_H2  = float(FH2_func(pH2_0,pCH4_0))
		## Soil-atmosphere exchanges
		QG    = max(Diff_z_X(T_surf,r_surf,a_eps,a_to,MG),0.)  / (max(z_vec)/len(z_vec)*1e4) * (pCH4_mol0-pCH4_mol_surf1) / (60*60*24)
		DiffG = -QG*rho0
		QH    = max(Diff_z_X(T_surf,r_surf,a_eps,a_to,MH),0.)  / (max(z_vec)/len(z_vec)*1e4) * (pH2_mol0-pH2_mol_surf1) / (60*60*24)
		DiffH = -QH*rho0
		## New state
		t_scaling = 80   # Can be adapted if the resulting dynamics are funky
		H2tot1  = H2tot0 + t_scaling*1/np.max([rho0,0.02])*(Volc - E*(pH2_0+pCH4_0/2) + DiffH + F_H2)*dt*(60*60*24)*S/Av
		CH4tot1 = max(CH4tot0 + t_scaling*1/np.max([rho0,0.02])*(DiffG + F_CH4)*dt*(60*60*24)*S/Av,1e-4*ntot0)

		ntot1   = ntot0 + (H2tot1-H2tot0) + (CH4tot1-CH4tot0)

		p1      = p0*ntot1/ntot0

		pH2_1  = H2tot1/ntot1
		pCH4_1 = CH4tot1/ntot1
		pCO2_1 = (1-(pCH4_1+pH2_1+pCO))*0.95
		pN2_1  = (1-(pCH4_1+pH2_1+pCO))*0.05

		T_av1       = T_func(p*1e-5,pH2_1,pCH4_1)
		rho1       = float(T_to_rho252(T_av1))    #Can be T_to_rho203, T_to_rho252, or T_to_rho273
		T_surf1    = max(T_av1,float(T_to_Tsurf252(T_av1)))    #Can be T_to_Tsurf203, T_to_Tsurf252, or T_to_Tsurf273
		T_vec1     = T_surf1 + a_T * z_vec

		Env1   = [pH2_1*p1*1e-5*alphaH*1e-3*Av,pCO2_1*p1*1e-5*np.exp(9345.17/T_surf1-167.8108+23.3585*np.log(T_surf1)+(0.023517-2.3656e-4*T_surf1+4.7036e-7*T_surf1**2)*35.0)*1e-3*Av,pCH4_1*p1*1e-5*alphaG*1e-3*Av,pN2_1*p1*1e-5*alphaN2*1e-3*Av]

		t_count += dt
		y1[np.where(y1 <= 1e-99)] = 1e-99
		y1[np.where(np.isnan(y1) == True)] = 1e-99

		if t_count > ratchet:
			t   += [t_count]
			ratchet = floor(t_count+ratchet_step)

			X0_t        = y1[(5*layers):(7*layers+layers)]
			NCT_t       = y1[(4*layers):(6*layers+layers)]
			HT_t        = y1[(0*layers):(0*layers+layers)]
			CT_t        = y1[(1*layers):(1*layers+layers)]
			NT_t        = y1[(3*layers):(5*layers+layers)]
			GT_t        = y1[(2*layers):(2*layers+layers)]

			NCT         = [NCT[i]       + [NCT_t[i]]        for i in np.arange(layers)]
			HT          = [HT[i]        + [HT_t[i]]         for i in np.arange(layers)]
			CT          = [CT[i]        + [CT_t[i]]         for i in np.arange(layers)]
			NT          = [NT[i]        + [NT_t[i]]         for i in np.arange(layers)]
			GT          = [GT[i]        + [GT_t[i]]         for i in np.arange(layers)]

			y		+= [y1.tolist()]
			pT      += [p1]
			T_surfT += [T_surf1]
			T_avT   += [T_av1]
			rhoT    += [rho1]
			pCH4T   += [pCH4_1]
			pH2T    += [pH2_1]
			pCO2T   += [pCO2_1]
			pN2T    += [pN2_1]
			DiffHT  += [DiffH]
			DiffGT  += [DiffG]

	t   += [t_count]
	ratchet = floor(t_count+ratchet_step)

	X0_t        = y1[(5*layers):(7*layers+layers)]
	NCT_t       = y1[(4*layers):(6*layers+layers)]
	HT_t        = y1[(0*layers):(0*layers+layers)]
	CT_t        = y1[(1*layers):(1*layers+layers)]
	NT_t        = y1[(3*layers):(5*layers+layers)]
	GT_t        = y1[(2*layers):(2*layers+layers)]

	NCT         = [NCT[i]       + [NCT_t[i]]        for i in np.arange(layers)]
	HT          = [HT[i]        + [HT_t[i]]         for i in np.arange(layers)]
	CT          = [CT[i]        + [CT_t[i]]         for i in np.arange(layers)]
	NT          = [NT[i]        + [NT_t[i]]         for i in np.arange(layers)]
	GT          = [GT[i]        + [GT_t[i]]         for i in np.arange(layers)]

	y		+= [y1.tolist()]
	pT      += [p1]
	T_surfT += [T_surf1]
	T_avT   += [T_av1]
	rhoT    += [rho1]
	pCH4T   += [pCH4_1]
	pH2T    += [pH2_1]
	pCO2T   += [pCO2_1]
	pN2T    += [pN2_1]
	DiffHT  += [DiffH]
	DiffGT  += [DiffG]

t_vec = np.append(np.arange(0,max(t),1000),max(t))
Slice = np.unique([np.argmin(abs(t-i)) for i in t_vec])

## Saves the output under the name that must be provided as argument
if len(Slice) <= 1:
	f = gzip.open(sys.argv[1]+'.results_dc','wb')
	pickle.dump([np.array(t),np.array(y),np.array(pT),np.array(T_surfT),np.array(T_avT),np.array(rhoT),np.array(pCH4T),np.array(pH2T),np.array(pCO2T),np.array(pN2T),np.array(DiffHT),np.array(DiffGT)],f)
	f.close()
else:
	f = gzip.open(sys.argv[1]+'.results_dc','wb')
	pickle.dump([np.array(t)[Slice],np.array(y)[Slice],np.array(pT)[Slice],np.array(T_surfT)[Slice],np.array(T_avT)[Slice],np.array(rhoT)[Slice],np.array(pCH4T)	[Slice],np.array(pH2T)[Slice],np.array(pCO2T)[Slice],np.array(pN2T)[Slice],np.array(DiffHT)[Slice],np.array(DiffGT)[Slice]],f)
	f.close()
