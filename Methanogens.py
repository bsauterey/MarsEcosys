############################################
###  THIS HERE DESCRIBES THE METABOLISM  ###
###  OF THE METHANOGENS                  ###
############################################

import numpy as np
from math import *

R = 8.314 
TS = 298


#### ENERGIES ####

# H,CO2,NH4,CH4,X,CO,CH3COOH,NO3,NO2,N2,H2SO4,H2S,XD
Catabolism = np.array([-1,-0.25,0,0.25,0,0,0,0,0,0,0,0,0])   # Catabolism stoichiometry
#Anabolism  = np.array([-2.1,-1,-0.2,0,1,0,0,0,0,0,0,0,0])    # Anabolism stoichiometry
Anabolism  = np.array([-2.4,-1,0,0,1,0,0,0,0,-0.1,0,0,0])    # Anabolism stoichiometry
    
NoC = 1                                      # Carbon source chain length
gamma = 4                                    # Carbon oxidation state in carbon source

# Standard Gibbs energy (J.mol eD-1)
deltaG0Cat = -32575                          
#deltaG0Cat = -36800                          
#deltaG0Ana = -12390 
deltaG0Ana = 28250 
#deltaG0Ana = -29200 

# Standard Enthalpy (J.mol eD-1)
deltaH0Cat = -63175
#deltaH0Cat = -58200 
#deltaH0Ana = -99700 
deltaH0Ana = 128000
#deltaH0Ana = -79900 #J/mol

def Dgdiss(NoC,gamma):
    """
    Returns the dissipaed energy during metabolism **REFERENCE**
    """
    return((200 + 18*(6 - NoC)**1.8 + exp(((-0.2 - gamma)**2)**0.16*(3.6 + 0.4*NoC)))*1000)

dgdiss = Dgdiss(NoC,gamma)

def DeltaG0(T,DeltaG0s,DeltaH0s):
    """
    Compute the DeltaG0 value when the temperature is different from the standard temperature TS = 298 K
    """
    global TS
    return(DeltaG0s*(T/TS)+DeltaH0s*((TS-T)/TS))

def DeltaG(T,DeltaG0s,C,S,DeltaH0s):
    """
    Computes the Gibbs energy of the reaction
    C is the list of reactants and products
    S is their stoichiometric coefficients
    (Must be in the same order : H,C,N,G,X!!)
    DeltaG0 and DeltaH0 are for T=298K
    """
    global R
    global TS
    
    DeltaG0T = DeltaG0(T,DeltaG0s,DeltaH0s)
    DeltaG = DeltaG0T + R*T*np.log(np.product(np.array(C)**np.array(S)))
    return(DeltaG)

def DeltaGcat(T,H,C,G):
    """
    Computes the catabolic gibbs energy for the H2 autotrophs
    """
    global deltaG0Cat
    global deltaH0Cat
    return(DeltaG(T,deltaG0Cat,[H,C,G],[-1,-0.25,0.25],deltaH0Cat))

def DeltaGana(T,H,C,N,x):
    """
    Computes the anabolic gibbs energy for the H2 autotrophs
    """
    global deltaG0Ana
    global deltaH0Ana

    return(DeltaG(T,deltaG0Ana,[H,C,N,x],[-2.4,-1,-0.1,1],deltaH0Ana))




#### RATES ####

def Yl(lam):
    """
    Returns the effective stoichiometry taking into account
    times catabolic reaction has to run to fuel anabolic reaction
    order is H C N G X
    """
    global Catabolism
    global Anabolism
    return(np.add(lam*Catabolism,Anabolism))


def Slim(vec,S): #TO OPTIMIZE TAKES TOO MUCH TIME
    """
    Finds the limiting substrate of a reaction
    """
    S2 = np.array(S)
    C2 = np.array(vec)
    stbalanced = np.abs(np.array(C2)/np.array(S2))
    return(vec[list(stbalanced).index(np.min(stbalanced))])

def Monod(qmax,k,x):
    """
    returns a reaction rate
    """
    return(qmax*(x/(k+x)))

def Mreq(mg,deltaGcat):
    """
    returns the maintenance energy requirement in term of catabolic energy
    """
    return(-mg/deltaGcat)

def QCat(deltaGcat,H,C,qmax,ks):
    """
    Computes the rate at which catabolic reaction occurs 
    """
    if deltaGcat < 0:
        return(Monod(qmax,ks,Slim([H,C],[-1,-0.25])))
    else:
        return(0)

def QMet(deltaGcat,qmax,ks,Smetlim):
    """
    Returns the rate at which the metabolism runs
    """
    if deltaGcat < 0:
        return(Monod(qmax,ks,Smetlim))
    else:
        return(0)

def QAna(deltaGcat,deltaGana,lam,qcat,qmet,mreq,qmax,ks,Sanalim):
    """
    Computes the rate at which anabolic reaction runs
    """
    global dgdiss

    if deltaGcat < 0 and lam > 0 and qcat > mreq:
        return((1/lam)*(qcat-mreq)*Monod(qmax,ks,Sanalim))/qmax
    else:
        return(0)

def Decay(mreq,qcat,deltaGcat,kd):
    """
    Computes the decay rate 
    """
    if mreq > qcat and deltaGcat <= 0:
        return(kd*(mreq-qcat)/mreq)
    elif deltaGcat > 0:
        return(kd)
    else:
        return(0)

def Gamma(thresh,slope,gmax,x,Qc):
    """
    Division function
    """
    if x > 2*thresh:
        gamma = 1/(1+np.exp(-slope*(np.log10((x-2*thresh)/thresh))))
    else:
        gamma = 0
    return(gamma)
