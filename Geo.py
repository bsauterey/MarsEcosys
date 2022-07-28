"""
This code is under MIT license. See the License.txt file.

This file contains the subroutines resolving the dynamics specific to the geoclimatic model alone.

Boris Sauterey
boris.sauterey@ens.fr
"""

from Constants import *
from math import *
from Technical_functions import *


def Run_Geo(Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,Cp,ALKo,ALKp,
            concCao,concCap,omegao0,omegap0,FHa,FCa,FGa,FCOa,FHoc,
            FCoc,FNoc,FGoc,FCOoc,FCH3COOHoc,FNO3oc,FNO2oc,FN2oc,FH2SO4oc,FH2Soc,FX0Doc,t,T0,T_0,fact,pCO2,Volc,tmax,
            t_exclu_m,FTnames,FTlist,t_intro_list,Delta,alphaC,Ccycle):

    Hatm0      = Hatm
    Catm0      = Catm
    Gatm0      = Gatm
    COatm0     = COatm
    Hoc0       = Hoc
    Coc0       = Coc
    Noc0       = Noc
    Goc0       = Goc
    COoc0      = COoc
    CH3COOHoc0 = CH3COOHoc
    NO3oc0     = NO3oc
    NO2oc0     = NO2oc
    N2oc0      = N2oc
    H2SO4oc0   = H2SO4oc
    H2Soc0     = H2Soc
    X0Doc0     = X0Doc
    T          = T_0
    Co0        = Co
    Cp0        = Cp
    ALKo0      = ALKo
    ALKp0      = ALKp
    concCao0   = concCao
    concCap0   = concCap

    t_temp     = 0
    delta_system = 0
    ti         = t
    
    # Choose the type of T° feedback (0: none; 1: atm. composition; 2: atm. composition and time):
    fT       = 1
    deltaT   = 0.2

    ## The 'while' ensures the geoclimatic conditions do not change to much
    ## so that the symplifying assumption of fixed biogenic fluxes remains not to false
    while (delta_system < Delta and t < tmax and (T - T_0) < deltaT):

        t0       = t

        a       = min(max(0.5,0.6 + 0.1*log10(Hatm/ntot/1e-4) - (log10(Gatm/ntot)+2)/15),1)
        Phot    = (2e10*(Catm/ntot/0.1)**(-0.2) * (Hatm/ntot/1e-4)**(-0.2) * (Gatm/ntot/1e-4)**a) * 365*60*60*24 * S /Av
        Conv    = Vco*COatm*p*1e-4 / (ntot*R*T)
        Esc     = -E*(Hatm+2*Gatm)/(ntot)
        CO2toCO = (1.8e10*(Catm/ntot/1e-1)**0.5 * (Hatm/ntot/1e-4)**0.4) * 365*60*60*24 * S /Av

        ## Compute the derivatives of the chemical elements assuming fix biogenic fluxes FX_is
        dHatm      = 2*Phot - CO2toCO + S * (FHa + Esc + Volc + Conv)
        dCatm      = - CO2toCO - Phot + S * (Conv + FCa)
        dGatm      = - Phot + S * (FGa)
        dCOatm     = CO2toCO + 2*Phot + S * (-Conv + FCOa)
        dHoc       = FHoc
        dNoc       = NO3oc*Fhyd*365 + NO2oc*Fhyd*365 + N2oc*Fhyd*365 + FNoc
#        dNoc       = FNoc
        dGoc       = FGoc
        dCOoc      = FCOoc
        dCH3COOHoc = FCH3COOHoc
        dNO3oc     = FNO3oc   - NO3oc*Fhyd*365
        dNO2oc     = FNO2oc   - NO2oc*Fhyd*365
        dN2oc      = FN2oc    - N2oc*Fhyd*365
        dH2SO4oc   = FH2SO4oc - H2SO4oc*Fhyd*365
        dH2Soc     = FH2Soc   - H2Soc*Fhyd*365
        dX0Doc     = FX0Doc


        if Ccycle == True: Catm,dCo,dCp,dALKo,dALKp,dconcCao,dconcCap,Coc1 = Carbon_Cycle(Catm,Co,Cp,ALKo,ALKp,concCao,concCap
                                                                                         ,omegao0,omegap0,T,T0,fact,
                                                                                          p,t,FCoc,alphaC)
        else: dCatm = 0; dCo = 0; dCp = 0; dALKo = 0; dALKp = 0; dconcCao = 0; dconcCap = 0; Coc1 = Coc


        if Hatm      <= 1e-20 and dHatm      < 0: dHatm      = 0 ; Hatm      = 1e-20
        if Catm      <= 1e-20 and dCatm      < 0: dCatm      = 0 ; Catm      = 1e-20    
        if Gatm      <= 1e-20 and dGatm      < 0: dGatm      = 0 ; Gatm      = 1e-20
        if COatm     <= 1e-20 and dCOatm     < 0: dCOatm     = 0 ; COatm     = 1e-20
        if Hoc       <= 1e-20 and dHoc       < 0: dHoc       = 0 ; Hoc       = 1e-20
        if Coc       <= 1e-20 and dCoc       < 0: dCoc       = 0 ; Coc       = 1e-20
        if Noc       <= 1e-20 and dNoc       < 0: dNoc       = 0 ; Noc       = 1e-20
        if Goc       <= 1e-20 and dGoc       < 0: dGoc       = 0 ; Goc       = 1e-20
        if COoc      <= 1e-20 and dCOoc      < 0: dCOoc      = 0 ; COoc      = 1e-20
        if CH3COOHoc <= 1e-20 and dCH3COOHoc < 0: dCH3COOHoc = 0 ; CH3COOHoc = 1e-20
        if NO3oc     <= 1e-20 and dNO3oc     < 0: dNO3oc     = 0 ; NO3oc     = 1e-20
        if NO2oc     <= 1e-20 and dNO2oc     < 0: dNO2oc     = 0 ; NO2oc     = 1e-20
        if N2oc      <= 1e-20 and dN2oc      < 0: dN2oc      = 0 ; N2oc      = 1e-20
        if H2SO4oc   <= 1e-20 and dH2SO4oc   < 0: dH2SO4oc   = 0 ; H2SO4oc   = 1e-20
        if H2Soc     <= 1e-20 and dH2Soc     < 0: dH2Soc     = 0 ; H2Soc     = 1e-20
        if Co        <= 1e-20 and dCo        < 0: dCo        = 0 ; Co        = 1e-20
        if Cp        <= 1e-20 and dCp        < 0: dCp        = 0 ; Cp        = 1e-20
        if ALKo      <= 1e-20 and dALKo      < 0: dALKo      = 0 ; ALKo      = 1e-20
        if ALKp      <= 1e-20 and dALKp      < 0: dALKp      = 0 ; ALKp      = 1e-20
        if concCao   <= 1e-20 and dconcCao   < 0: dconcCao   = 0 ; concCao   = 1e-20
        if concCap   <= 1e-20 and dconcCap   < 0: dconcCap   = 0 ; concCap   = 1e-20

        ## Adjust the geoclimatic timestep
        vec        = np.array([Hatm,Co,Cp,ALKo,ALKp,concCao,concCap,Catm,Gatm,COatm,Hoc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc
                               ,N2oc,H2SO4oc,H2Soc,X0Doc])
        der_vec    = np.array([dHatm,dCo,dCp,dALKo,dALKp,dconcCao,dconcCap,dCatm,dGatm,dCOatm,dHoc,dNoc,dGoc,dCOoc,dCH3COOHoc
                               ,dNO3oc,dNO2oc,dN2oc,dH2SO4oc,dH2Soc,dX0Doc])
        dt         = min([min(abs(vec[np.where(der_vec != 0)[0]]
                              /der_vec[np.where(der_vec != 0)[0]]*Delta/10)),tmax-t,t_exclu_m-t,deltaT*1e8]+
                         (np.array([t_intro_list[i] for i in FTnames])-t).tolist())
        
        var1 = np.array([Hatm,Catm,Gatm,COatm,Hoc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,
                Cp,ALKo,ALKp,concCao,concCap])
        dvar = np.array([dHatm,dCatm,dGatm,dCOatm,dHoc,dNoc,dGoc,dCOoc,dCH3COOHoc,dNO3oc,dNO2oc,dN2oc,dH2SO4oc,dH2Soc,dX0Doc,dCo,
                dCp,dALKo,dALKp,dconcCao,dconcCap])
        var2 = (var1 * np.exp(dvar/var1*dt)).tolist()
        
        [Hatm1,Catm1,Gatm1,COatm1,Hoc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1,Co1,        
         Cp1,ALKo1,ALKp1,concCao1,concCap1] = var2 

        #print(var2)
        #print(Catm1,dCatm,dt)

        t        = t + dt
        t1       = t
        t_temp   = t_temp + dt

        if   fT == 0: T = T
        elif fT == 1: T = T0 + 19*np.log10(Catm/ntot*1e6/pCO2) + 5*np.log10(1+Gatm1/ntot*1e6/3)
        #elif fT == 1: T = 253.89 + 77.67*np.sqrt(Catm/ntot) + 5*np.log10(1+Gatm1/ntot*1e6/3)
        elif fT == 2: T = T0 + 19*np.log10(Catm/ntot*1e6/pCO2) + 5*np.log10(1+Gatm1/ntot*1e6/3) + 10*(t*1e-9)

        a       = min(max(0.5,0.6 + 0.1*log10(Hatm1/ntot/1e-4) - (log10(Gatm1/ntot)+2)/15),1)
        Phot    = (2e10*(Catm1/ntot/0.1)**(-0.2) * (Hatm1/ntot/1e-4)**(-0.2) * (Gatm1/ntot/1e-4)**a) * 365*60*60*24 * S /Av
        Conv    = Vco*COatm1*p*1e-4 / (ntot*R*T)
        Esc     = -E*(Hatm1+2*Gatm1)/(ntot)
        CO2toCO = (1.8e10*(Catm1/ntot/1e-1)**0.5 * (Hatm1/ntot/1e-4)**0.4) * 365*60*60*24 * S /Av

        ## Compute the derivatives of the chemical elements assuming fix biogenic fluxes FX_is
        dHatm      = 2*Phot - CO2toCO + S * (FHa + Esc + Volc + Conv)
        dCatm      = - CO2toCO - Phot + S * (Conv + FCa)
        dGatm      = - Phot + S * (FGa)
        dCOatm     = CO2toCO + 2*Phot + S * (-Conv + FCOa)
        dHoc       = FHoc
        dNoc       = NO3oc*Fhyd*365 + NO2oc*Fhyd*365 + N2oc*Fhyd*365 + FNoc
#        dNoc       = FNoc
        dGoc       = FGoc
        dCOoc      = FCOoc
        dCH3COOHoc = FCH3COOHoc
        dNO3oc     = FNO3oc   - NO3oc1*Fhyd*365
        dNO2oc     = FNO2oc   - NO2oc1*Fhyd*365
        dN2oc      = FN2oc    - N2oc1*Fhyd*365
        dH2SO4oc   = FH2SO4oc - H2SO4oc1*Fhyd*365
        dH2Soc     = FH2Soc   - H2Soc1*Fhyd*365
        dX0Doc     = FX0Doc

        if Ccycle == True: Catm1,dCo,dCp,dALKo,dALKp,dconcCao,dconcCap,Coc2 = Carbon_Cycle(Catm1,Co1,Cp1,ALKo1,ALKp1,concCao1
                                                                    ,concCap1,omegao0,omegap0,T,T0,fact,p,t,FCoc,alphaC)
        else: dCatm = 0; dCo = 0; dCp = 0; dALKo = 0; dALKp = 0; dconcCao = 0; dconcCap = 0; Coc2 = Coc1

        ## Adjust the geoclimatic timestep
        vec     = np.array([Hatm1,Co1,Cp1,ALKo1,ALKp1,concCao1,concCap1,Catm1,Gatm1,COatm1,Hoc1,Noc1,Goc1,COoc1,CH3COOHoc1,
                            NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1])
        der_vec = np.array([dHatm,dCo,dCp,dALKo,dALKp,dconcCao,dconcCap,dCatm,dGatm,dCOatm,dHoc,dNoc,dGoc,dCOoc,dCH3COOHoc,
                            dNO3oc,dNO2oc,dN2oc,dH2SO4oc,dH2Soc,dX0Doc])
        dt         = min([min(abs(vec[np.where(der_vec != 0)[0]]
                              /der_vec[np.where(der_vec != 0)[0]]*Delta/10)),tmax-t,t_exclu_m-t,deltaT*1e8]+
                         (np.array([t_intro_list[i] for i in FTnames])-t).tolist())
        
        var1 = np.array([Hatm1,Catm1,Gatm1,COatm1,Hoc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1,Co1,
                Cp1,ALKo1,ALKp1,concCao1,concCap1])
        dvar = np.array([dHatm,dCatm,dGatm,dCOatm,dHoc,dNoc,dGoc,dCOoc,dCH3COOHoc,dNO3oc,dNO2oc,dN2oc,dH2SO4oc,dH2Soc,dX0Doc,dCo,
                dCp,dALKo,dALKp,dconcCao,dconcCap])
        var2 = (var1 * np.exp(dvar/var1*dt)).tolist()
        [Hatm2,Catm2,Gatm2,COatm2,Hoc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2,Co2,        
         Cp2,ALKo2,ALKp2,concCao2,concCap2] = var2    
        
        t       = t + dt
        t2      = t
        t_temp  = t_temp + dt

        if   fT == 0: T = T
        elif fT == 1: T = T0 + 19*np.log10(Catm/ntot*1e6/pCO2) + 5*np.log10(1+Gatm1/ntot*1e6/3)
        #elif fT == 1: T = 253.89 + 77.67*np.sqrt(Catm/ntot) + 5*np.log10(1+Gatm1/ntot*1e6/3)
        elif fT == 2: T = T0 + 19*np.log10(Catm/ntot*1e6/pCO2) + 5*np.log10(1+Gatm1/ntot*1e6/3) + 10*(t*1e-9)

        [Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,Cp,ALKo
         ,ALKp,concCao,concCap] = oscillation_geo(Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,
                                                  H2SO4oc,H2Soc,X0Doc,Co,Cp,ALKo,ALKp,concCao,concCap,Hatm1,Catm1,Gatm1,
                                                  COatm1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,
                                                  H2Soc1,X0Doc1,Co1,Cp1,ALKo1,ALKp1,concCao1,concCap1,Hatm1,Catm1,Gatm1,COatm1,
                                                  Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,
                                                  X0Doc1,Co1,Cp1,ALKo1,ALKp1,concCao1,concCap1,t0,t1,t2)

        if   fT == 0: T = T
        elif fT == 1: T = T0 + 19*np.log10(Catm/ntot*1e6/pCO2) + 5*np.log10(1+Gatm1/ntot*1e6/3)
        #elif fT == 1: T = 253.89 + 77.67*np.sqrt(Catm/ntot) + 5*np.log10(1+Gatm1/ntot*1e6/3)
        elif fT == 2: T = T0 + 19*np.log10(Catm/ntot*1e6/pCO2) + 5*np.log10(1+Gatm1/ntot*1e6/3) + 10*(t*1e-9)

        var  = [Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc]
        var0 = [Hatm0,Catm0,Gatm0,COatm0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0]
        dvar = abs(1 - np.array(var)/np.array(var0))

        delta_system = max(dvar)

        for i in FTnames:
            if FTlist[i] is False:
                if t >= t_intro_list[i]:
                    delta_system = Delta*2
        
    return(Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,
           Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0,t,T)



def Run_abio(Hatm,tmax=1e4,Volc=Volc):

    # The goal of this is to obtain the prebiotic partial pressure of H2
    dHatm  = Hatm
    HatmT  = [Hatm]
    CatmT  = [Catm]
    GatmT  = [Gatm]
    COatmT = [COatm]
    t      = [0]
    #while dHatm > 1e-10*HatmT[-1]:
    while t[-1] < tmax:
        
        ## Phototochemical reactions
        a       = min(max(0.5,0.6 + 0.1*log10(HatmT[-1]/ntot/1e-4) - (log10(GatmT[-1]/ntot)+2)/15),1)
        Phot    = (2e10*(Catm/ntot/0.1)**(-0.2) * (HatmT[-1]/ntot/1e-4)**(-0.2) * (GatmT[-1]/ntot/1e-4)**a) * 365*60*60*24 * S /Av
        Conv    = Vco*COatmT[-1]*p*1e-4 / (ntot*R*T)
        Esc     = -E*(HatmT[-1]+2*GatmT[-1])/(ntot)
        CO2toCO = (1.8e10*(Catm/ntot/1e-1)**0.5 * (HatmT[-1]/ntot/1e-4)**0.4) * 365*60*60*24 * S /Av

        ## Compute the derivatives of the chemical elements assuming fix biogenic fluxes FX_is
        dHatm      = 2*Phot - CO2toCO + S * (Esc + Volc + Conv)
        dCatm      = - CO2toCO - Phot + S * Conv
        dGatm      = - Phot
        dCOatm     = CO2toCO + 2*Phot + S * (-Conv)

        vec        = np.array([Hatm,Gatm,COatm])
        der_vec    = np.array([dHatm,dGatm,dCOatm])
        dt         = 10

        HatmT.append(max(1e-20*1e-6*ntot,HatmT[-1]  * np.exp(dHatm/HatmT[-1]*dt)))
        CatmT.append(max(1e-20*1e-6*ntot,CatmT[-1]  * np.exp(dCatm/CatmT[-1]*dt)))
        GatmT.append(max(1e-20*1e-6*ntot,GatmT[-1]  * np.exp(dGatm/GatmT[-1]*dt)))
        COatmT.append(max(1e-20*1e-6*ntot,COatmT[-1]  * np.exp(dCOatm/COatmT[-1]*dt)))

        t.append(t[-1] + dt)

    HatmT  = np.array(HatmT)
    CatmT  = np.array(CatmT)
    GatmT  = np.array(GatmT)
    COatmT = np.array(COatmT)
    t      = np.array(t)

    return(t,HatmT,GatmT,COatmT)



def initialize_carbon_cycle(Catm,T,T0):

    ###compute Omegao0 and Omegap0
    K1o       = 10**(-17.788+0.073104*T0+0.0051087*35-1.1463e-4*T0**2)
    K2o       = 10**(-20.919+0.064209*T0+0.011887*35-8.7313e-5*T0**2)
    HCO2      = 1e-5*exp(9345.17/T0-167.8108+23.3585*log(T0)+(0.023517-2.3656e-4*T0+4.7036e-7*T0**2)*35)
    Kso       = 10**(-171.9065-0.077993*T0+2839.319/T0+71.595*log10(T0)
                     +(-0.77712+0.0028426*T0+178.34/T0)*35**0.5-0.0711*35.0+0.0041249*35**1.5)

    Tdeep     = max(min(1.1*T0-40.0,T0),271.15)
    Sthick    = 700*(1-(1-fsed)*tini/4)
    Qfactor   = (1-(3.8-tini)/4.5)**(-nout)
    Tpore     = Tdeep+Qfactor*Sthick/K
    K1p       = 10**(-17.788+0.073104*Tpore+0.0051087*35-1.1463e-4*Tpore**2)
    K2p       = 10**(-20.919+0.064209*Tpore+0.011887*35-8.7313e-5*Tpore**2)
    Ksp       = 10**(-171.9065-0.077993*Tpore+2839.319/Tpore+71.595*log10(Tpore)
                     +(-0.77712+0.0028426*Tpore+178.34/Tpore)*35**0.5-0.0711*35+0.0041249*35**1.5)

    concCO2o  = 280e-6*1e5*HCO2
    concCO2p  = concCO2o
    concHCO3o = concCO2o*K1o/concHo0
    concCO3o  = concHCO3o*K2o/concHo0
    Co        = (concCO2o+concHCO3o+concCO3o+so*280e-6*p)*Mo
    ALKo      = (2*concCO3o+concHCO3o)*Mo
    concCao   = ALKo/Mo+(1e-14/concHo0-concHo0)
    omegao0   = concCao*concCO3o/Kso

    Cp        = Co
    ALKp      = ALKo
    coef1     = ALKp/K1p/K2p/Mp
    coef2     = (ALKp-Cp)/K2p/Mp
    coef3     = (ALKp-2*Cp)/Mp
    delta     = coef2**2-4*coef1*coef3
    delta     = max(0.0,delta)
    concHp0   = (-coef2+delta**0.5)/(2*coef1)
    concHp0   = max(concHp0,1e-14)
    concHCO3p = concCO2p*K1p/concHp0
    concCO3p  = concHCO3p*K2p/concHp0
    concCap   = ALKp/Mp+(1e-14/concHp0-concHp0)
    omegap0   = concCap*concCO3p/Ksp

    ###compute Co and Cp
    K1o       = 10**(-17.788+0.073104*T+0.0051087*35-1.1463e-4*T**2)
    K2o       = 10**(-20.919+0.064209*T+0.011887*35-8.7313e-5*T**2)
    HCO2      = 1e-5*exp(9345.17/T-167.8108+23.3585*log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35)
    Kso       = 10**(-171.9065-0.077993*T+2839.319/T+71.595*log10(T)
                     +(-0.77712+0.0028426*T+178.34/T)*35**0.5-0.0711*35.0+0.0041249*35**1.5)
    Tdeep     = max(min(1.1*T-40.0,T),271.15)
    Sthick    = 700
    Qfactor   = (1-3.8/4.5)**(-nout)
    Tpore     = Tdeep+Qfactor*Sthick/K
    K1p       = 10**(-17.788+0.073104*Tpore+0.0051087*35-1.1463e-4*Tpore**2)
    K2p       = 10**(-20.919+0.064209*Tpore+0.011887*35-8.7313e-5*Tpore**2)
    Ksp       = 10**(-171.9065-0.077993*Tpore+2839.319/Tpore+71.595*log10(Tpore)
                     +(-0.77712+0.0028426*Tpore+178.34/Tpore)*35**0.5-0.0711*35+0.0041249*35**1.5)

    concHo    = concHo0_hadean
    concCO2o  = Catm/ntot*p*HCO2
    concCO2p  = concCO2o
    counter   = 0

    while counter < 100:
        concHCO3o = concCO2o*K1o/concHo
        concCO3o  = concHCO3o*K2o/concHo
        Co        = Mo*(concCO2o+concHCO3o+concCO3o)+Catm
        ALKo      = Mo*(2*concCO3o+concHCO3o)
        coef1     = ALKo/K1o/K2o*(1+so/HCO2)/Mo
        coef2     = (ALKo-Co)/K2o/Mo
        coef3     = (ALKo-2*Co)/Mo
        delta     = coef2**2-4*coef1*coef3
        delta     = max(0.0,delta)
        concHo    = (-coef2+delta**0.5)/(2*coef1)
        counter   = counter+1

    concCao   = 1.0*Kso/concCO3o

    Cp        = Co
    ALKp      = ALKo
    coef1     = ALKp/K1p/K2p/Mp
    coef2     = (ALKp-Cp)/K2p/Mp
    coef3     = (ALKp-2*Cp)/Mp
    delta     = coef2**2-4*coef1*coef3
    delta     = max(0.0,delta)
    concHp0   = (-coef2+delta**0.5)/(2*coef1)
    concHp0   = max(concHp0,1e-14)
    concHCO3p = concCO2p*K1p/concHp0
    concCO3p  = concHCO3p*K2p/concHp0
    concCap   = 1.0*Ksp/concCO3p

    return(Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0)



def Carbon_Cycle(Catm,Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0,T,T0,fact,p,time,FCoc,alphaC):
    
    T1 = T0
    #tga=(tini-time*1e-9) 	#time in Ga
    tga=tini
    Qfactor=(1-(3.8-tga)/4.5)**(-nout)
    rspreading=Qfactor**beta
    Sthick=700*(1-(1-fsed)*tga/4)
    Tdeep=max(min(1.1*T-40,T),271.15)
    Tpore=Tdeep+Qfactor*Sthick/K
    fland=max(0, 1-1/(1/(1-Larchean)+exp(-10*(tga-tgrowth))))
    fbio=max(0, 1-1/(1/(1-Bprecambrian)+exp(-10*(tga-tbio))))
    #fbio=1
    Kso=10**(-171.9065-0.077993*T+2839.319/T+71.595*log10(T)
             +(-0.77712+0.0028426*T+178.34/T)*35**0.5-0.0711*35.0+0.0041249*35**1.5)
    Ksp=10**(-171.9065-0.077993*Tpore+2839.319/Tpore+71.595*log10(Tpore)
             +(-0.77712+0.0028426*Tpore+178.34/Tpore)*35**0.5-0.0711*35.0+0.0041249*35**1.5)
    
    concHp,concCO2p,concHCO3p,concCO3p,none=ocean_chemistry(Cp,ALKp,Mp,Tpore,p,0.0,alphaC)
    concHo,concCO2o,concHCO3o,concCO3o,Catmnew=ocean_chemistry(Co,ALKo,Mo,T,p,so,alphaC)

#    concCao= 0.5*(ALKo-ALKoini)+concCaoini
#    concCap= 0.5*(ALKp-ALKpini)+concCapini

    omegao=concCao*concCO3o/Kso
    omegap=concCap*concCO3p/Ksp

    Fout=0
    Fcarb=0
    Fsil=0
    Fdiss=0
    Pocean=0
    Ppore=0
    #Fout=Fout0*Qfactor**mout
    Fout=(Fout0*Qfactor**mout)/fact
    Fcarb=fland*Fcarb0*(Catmnew/ntot/280e-6)**zeta*exp((T-T0)/Te)
    Fsil=fbio*fland*Fsil0*(Catmnew/ntot/280e-6)**zeta*exp((T-T0)/Te) # The Geo-To-Bio delay seems to be linked to this
    Fdiss=Fdiss0*rspreading*exp(-Ebas/R*(1/Tpore-1/Tpore0))*(concHp/concHo0)**gam_ma
    Pocean=Pocean0*fland*max(0,omegao-1)**n/max(1e-30,abs(omegao0-1))**n
    Ppore=Ppore0*max(0,omegap-1)**n/max(1e-30,abs(omegap0-1))**n
    dCo   = -J*(Co-Cp)/Mo + (Fout + Fcarb - Pocean) + FCoc*(DOZ/10*S)
    dCp   =  J*(Co-Cp)/Mp - Ppore
    dALKo = -J*(ALKo-ALKp)/Mo + (2*Fsil + 2*Fcarb - 2*Pocean)
    dALKp = J*(ALKo-ALKp)/Mp + (2*Fdiss - 2*Ppore)
    dconcCao= 0.5*dALKo/Mo
    dconcCap= 0.5*dALKp/Mp

    return(Catmnew,dCo,dCp,dALKo,dALKp,dconcCao,dconcCap,concCO2o)



def ocean_chemistry(C,ALK,M,T,p,s1,alphaC):

    #T = T0
    K1=10**(-17.788+0.073104*T+0.0051087*35-1.1463e-4*T**2)
    K2=10**(-20.919+0.064209*T+0.011887*35-8.7313e-5*T**2)
    K1=10**(-17.788+0.073104*273+0.0051087*35-1.1463e-4*273**2)
    K2=10**(-20.919+0.064209*273+0.011887*35-8.7313e-5*273**2)
    HCO2=1e-5*exp(9345.17/T-167.8108+23.3585*log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35.0)
    HCO2=1e-5*alphaC

    coef1=ALK/K1/K2*(1+s1/HCO2)/M
    coef2=(ALK-C)/K2/M
    coef3=(ALK-2*C)/M
    delta=coef2**2-4*coef1*coef3
    delta=max(0.0, delta)

    concH=(-coef2+delta**0.5)/(2*coef1)
    concH=max(concH,1e-14)
    Catm=C/(K1*K2*HCO2/concH**2+K1*HCO2/concH+HCO2+s1)*s1
    concCO2=Catm/ntot*p*HCO2
    concHCO3=concCO2*K1/concH
    concCO3=concHCO3*K2/concH

    return(concH,concCO2,concHCO3,concCO3,Catm)