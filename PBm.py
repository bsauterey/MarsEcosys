"""
This code is under MIT license. See the License.txt file.

This file is the main routine of the Planet-Biosphere model (PBm) that resolves numerically the dynamical coupling between a planet (ocean, atmosphere, core) and its biosphere.

Boris Sauterey
boris.sauterey@ens.fr
"""

## Python libraries
import numpy as np
import random
import time as libtime
import pickle
import scipy as sc
import matplotlib.pyplot as plt
from IPython import display
from scipy import stats
from Constants import *
from math import *

## Home-made libraries
import Methanogens as mt
import NO3Methanotrophs as NO3mtt
import NO2Methanotrophs as NO2mtt
import H2SO4Methanotrophs as H2SO4mtt
import Acetogens1 as ac
import Acetogens2 as ac2
import Acetotrophs as at
import Fermentors as f
import PhotoH2 as pH2
from Geo import *
from Bio import *
from Technical_functions import *

def Run_coupled_model(tmax,FTnames,FTlist,t_intro_list,T0 = T,pCO2 = pCO2,Volc = Volc,t_exclu_m = 1e100,Input = [],
                      Ccycle = False, extinction = True, biotmax = 365, fact = 1):

    global Catm
    global NO3inf
    global NO2inf
    global H2SO4inf
    global alphaC
    global T
    dum = True

#    print(t_intro_list)
    # Prebiotic atmosphere
    Hatm           = ntot*200*1e-6
    t,HatmT,GatmT,COatmT = Run_abio(Hatm,1e7)
    Catm           = pCO2 * ntot * 1e-6
    Hatm           = HatmT[-1]
    Gatm           = 1e-14  * ntot * 1e-6
    COatm          = COatmT[-1]
    # Prebiotic surface ocean
    Hinf           = Hatm/ntot  * alphaH
    Cinf           = Catm/ntot  * alphaC
    Ninf           = 1e-7
    Ginf           = Gatm/ntot  * alphaG
    COinf          = COatm/ntot * alphaCO
    CH3COOHinf     = 1e-20
    N2inf          = N2atm/ntot * alphaN2
    X0Dinf         = 1e-20
    H2Sinf         = 1e-20

    # Initial temperature
    T              = T0
    # Prebiotic deep ocean
    Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0=initialize_carbon_cycle(Catm,T,T0)

    #print(T)
    #print(Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0)

    # Prepare all the vectors needed later
    ## Fluxes
    FHa_vec         = [0]
    FCa_vec         = [0]
    FGa_vec         = [0]
    FCOa_vec        = [0]
    FHoc_vec        = [0]
    FCoc_vec        = [0]
    FNoc_vec        = [0]
    FGoc_vec        = [0]
    FCOoc_vec       = [0]
    FN2oc_vec       = [0]
    FNO3oc_vec      = [0]
    FNO2oc_vec      = [0]
    FH2SO4oc_vec    = [0]
    FH2Soc_vec      = [0]
    FCH3COOHoc_vec  = [0]
    FX0Doc_vec      = [0]
    NPP_vec         = [0]
    NPP_vec_list    = {'meth':[],'NO3metht':[],'NO2metht':[],'H2SO4metht':[],'acet':[],'acet2': [],'acett':[],'ferm':[],
                       'photoH2':[]}

    ## Chemical species
    ### Atm:
    HatmT          = [Hatm]
    CatmT          = [Catm]
    GatmT          = [Gatm]
    COatmT         = [COatm]
    N2atmT         = [N2atm]
    ### Surface ocean:
    N_vec          = [Ninf]
    CH3COOH_vec    = [CH3COOHinf]
    C_vec          = [Cinf]
    Co_vec         = [Co]
    NO3_vec        = [NO3inf]
    NO2_vec        = [NO2inf]
    H2SO4_vec      = [H2SO4inf]
    H2S_vec        = [H2Sinf]
    CO_vec         = [COinf]
    N2_vec         = [N2inf]
    X0D_vec        = [X0Dinf]
    ### Deep ocean:
    HocT           = [Hinf]
    CocT           = [Cinf]
    NocT           = [Ninf]
    GocT           = [Ginf]
    COocT          = [COinf]
    CH3COOHocT     = [CH3COOHinf]
    NO3ocT         = [NO3inf]
    NO2ocT         = [NO2inf]
    N2ocT          = [N2inf]
    H2SO4ocT       = [H2SO4inf]
    H2SocT         = [H2Sinf]
    X0DocT         = [X0Dinf]
    CoT            = [Co]
    CpT            = [Cp]
    ALKoT          = [ALKo]
    ALKpT          = [ALKp]
    concCaoT       = [concCao]
    concCapT       = [concCap]
    ## Other stuff
    time_vec        = [0]
    T_vec           = [T]
    time_atm        = [0]

    FHa            = 0
    FGa            = 0
    FCOa           = 0
    FN2a           = 0
    FHoc           = 0
    FCoc           = 0
    FNoc           = 0
    FGoc           = 0
    FCOoc          = 0
    FCH3COOHoc     = 0
    FNO3oc         = 0
    FNO2oc         = 0
    FH2SO4oc       = 0
    FH2Soc         = 0
    FN2oc          = 0
    FX0Doc         = 0

    init_FT_list,FTlist,t_intro_list = init_bio(FTnames,FTlist,t_intro_list,T)
    Biomass_list     = {'meth':[[],[]],'NO3metht':[[],[]],'NO2metht':[[],[]],'H2SO4metht':[[],[]],'acet':[[],[]],
                        'acet2': [[],[]],'acett':[[],[]],'ferm':[[],[]],'photoH2':[[],[]]}

    clock_intro_list = {'meth':0,'NO3metht':0,'NO2metht':0,'H2SO4metht':0,'acet':0,'acet2':0,'acett':0,'ferm':0,'photoH2':0}
    NPP_list         = {'meth':0,'NO3metht':0,'NO2metht':0,'H2SO4metht':0,'acet':0,'acet2':0,'acett':0,'ferm':0,'photoH2':0}


    t_0 = libtime.time()

    index_slope       = 0

    count      = 0
    count_dump = 1
    time0      = 0
    Catm0      = CatmT[-1]
    Hatm0      = HatmT[-1]
    Gatm0      = GatmT[-1]
    COatm0     = COatmT[-1]
    T_0        = T_vec[-1]
    Hoc0       = HocT[-1]
    Coc0       = CocT[-1]
    Noc0       = NocT[-1]
    Goc0       = GocT[-1]
    COoc0      = COocT[-1]
    CH3COOHoc0 = CH3COOHocT[-1]
    NO3oc0     = NO3ocT[-1]
    NO2oc0     = NO2ocT[-1]
    N2oc0      = N2ocT[-1]
    H2SO4oc0   = H2SO4ocT[-1]
    H2Soc0     = H2SocT[-1]
    X0Doc0     = X0DocT[-1]

    Hinf       = Hatm0/ntot  * alphaH
    Ginf       = Gatm0/ntot  * alphaG
    COinf      = COatm0/ntot * alphaCO

    while time0 < tmax:

        t0 = libtime.time()
        Delta = 0.1
        flat = False
        if len(HatmT)-index_slope > 10:
            flat,Delta = Delta_system(HatmT,GatmT,COatmT,HocT,CocT,NocT,GocT,COocT,CH3COOHocT,NO3ocT,NO2ocT,N2ocT,
                                      H2SO4ocT,H2SocT,X0DocT,time_atm)

        if sum(FTlist[i] for i in FTnames) == 0:
            Delta = 0.1

        ################################################
        # FIRST SUBSTEP:
        ############

        Env        = np.array([Hinf,Cinf,Ninf,Ginf,COinf,CH3COOHinf,NO3inf,NO2inf,N2inf,H2SO4inf,H2Sinf,X0Dinf,T,
                      Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0])
        if time_atm[-1] == 0: init_chem  = Env

        ####################################
        # Bio-module => Compute the bioresulting geoclimatic derivatives
        ############

        HT,CT,NT,GT,COT,CH3COOHT,NO3T,NO2T,N2T,H2SO4T,H2ST,X0DT,Bio,NPPT,NPPT_list,dX0DT,t,FTlist,Recycling = \
        Run_local_biology(Env,init_chem,FTnames,init_FT_list,biotmax,FTlist,time0,extinction=extinction)

        init_chem  = np.array([HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],
                              NO3T[-1],NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1],X0DT[-1]])

        for i in FTnames:
            if FTlist[i] is True:
                init_FT_list[i] = np.array([Bio[i][0][-1],Bio[i][1][-1]])

        index = min(min(np.where(t > max(t)*0.5)))
        # Global biogenic fluxes:
        FHa         = -QH  * (Hinf-sum(HT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10
        FCa         = -QC  * (Cinf-sum(CT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10
        FGa         = -QG  * (Ginf-sum(GT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10
        FCOa        = -QCO * (COinf-sum(COT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10

        FHoc        = -Fconv * (Hoc0-sum(HT[index:]            *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FCoc        = -Fconv * (Coc0-sum(CT[index:]            *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FNoc        = -Fconv * (Noc0-sum(NT[index:]            *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FGoc        = -Fconv * (Goc0-sum(GT[index:]            *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FCOoc       = -Fconv * (COoc0-sum(COT[index:]          *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FCH3COOHoc  = -Fconv * (CH3COOHoc0-sum(CH3COOHT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FNO3oc      = -Fconv * (NO3oc0-sum(NO3T[index:]        *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FNO2oc      = -Fconv * (NO2oc0-sum(NO2T[index:]        *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FN2oc       = -Fconv * (N2oc0-sum(N2T[index:]          *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FH2SO4oc    = -Fconv * (H2SO4oc0-sum(H2SO4T[index:]    *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FH2Soc      = -Fconv * (H2Soc0-sum(H2ST[index:]        *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FX0Doc      = -Fconv * (X0Doc0-sum(X0DT[index:]        *diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ

        NPP         = sum(NPPT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:]))*365*10*S
        for i in FTnames:
            if FTlist[i] is True:
                NPP_list[i] = sum(NPPT_list[i][index:]*diff(t[index-1:]))/sum(diff(t[index-1:]))*365*10*S

        ####################################
        # Geo-module
        ############
        Hatm1,Catm1,Gatm1,COatm1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1, \
        Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0,time1,T1 = \
        Run_Geo(Hatm0,Catm0,Gatm0,COatm0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0,
                Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0,FHa,FCa,FGa,FCOa,FHoc,FCoc,FNoc,FGoc,FCOoc,FCH3COOHoc,
                FNO3oc,FNO2oc,FN2oc,FH2SO4oc,FH2Soc,FX0Doc,time0,T0,T,fact,pCO2,Volc,tmax,
                t_exclu_m,FTnames,FTlist,t_intro_list,Delta,alphaC,Ccycle)

        alphaC = exp(9345.17/T1-167.8108+23.3585*log(T1)+(0.023517-2.3656e-4*T1+4.7036e-7*T1**2)*35.0)

        t1 = libtime.time()

        # Updates the local oceanic environment for the biological species
        Hinf       = Hatm1/ntot  * alphaH
        Cinf       = Catm1/ntot  * alphaC
        Ninf       = Noc1
        Ginf       = Gatm1/ntot  * alphaG
        COinf      = COatm1/ntot * alphaCO
        N2inf      = N2atm/ntot  * alphaN2
        Env        = np.array([Hinf,Cinf,Ninf,Ginf,COinf,CH3COOHinf,NO3inf,NO2inf,N2inf,H2SO4inf,H2Sinf,X0Dinf,T1,
                      Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1])

        ################################################
        # SECOND SUBSTEP:
        ############

        ####################################
        # Bio-module => Compute the bioresulting geoclimatic derivatives
        ############

        HT,CT,NT,GT,COT,CH3COOHT,NO3T,NO2T,N2T,H2SO4T,H2ST,X0DT,Bio,NPPT,NPPT_list,dX0DT,t,FTlist,Recycling = \
        Run_local_biology(Env,init_chem,FTnames,init_FT_list,biotmax,FTlist,time1,extinction=extinction)

        init_chem  = np.array([HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],
                               NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1],X0DT[-1]])

        index = min(min(np.where(t > max(t)*0.5)))

        # Global biogenic fluxes:
        FHa         = -QH    * (Hinf-sum(HT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10
        FCa         = -QG    * (Cinf-sum(CT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10
        FGa         = -QG    * (Ginf-sum(GT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10
        FCOa        = -QCO   * (COinf-sum(COT[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:])))*365*10

        FHoc        = -Fconv * (Hoc1-sum(HT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FCoc        = -Fconv * (Coc1-sum(CT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FNoc        = -Fconv * (Noc1-sum(NT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FGoc        = -Fconv * (Goc1-sum(GT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FCOoc       = -Fconv * (COoc1-sum(COT[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FCH3COOHoc  = -Fconv * (CH3COOHoc1-sum(CH3COOHT[index:] * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FNO3oc      = -Fconv * (NO3oc1-sum(NO3T[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FNO2oc      = -Fconv * (NO2oc1-sum(NO2T[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FN2oc       = -Fconv * (N2oc1-sum(N2T[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FH2SO4oc    = -Fconv * (H2SO4oc1-sum(H2SO4T[index:]     * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FH2Soc      = -Fconv * (H2Soc1-sum(H2ST[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ
        FX0Doc      = -Fconv * (X0Doc1-sum(X0DT[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:]))) / DOZ

        NPP         = sum(NPPT[index:]*diff(t[index-1:]))/sum(diff(t[index-1:]))*365*10*S
        for i in FTnames:
            if FTlist[i] is True:
                NPP_list[i] = sum(NPPT_list[i][index:]*diff(t[index-1:]))/sum(diff(t[index-1:]))*365*10*S

        ####################################
        # Geo-module
        ############

        Hatm2,Catm2,Gatm2,COatm2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2, \
        Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0,time2,T2 = \
        Run_Geo(Hatm1,Catm1,Gatm1,COatm1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1,
                Co,Cp,ALKo,ALKp,concCao,concCap,omegao0,omegap0,FHa,FCa,FGa,FCOa,FHoc,FCoc,FNoc,FGoc,FCOoc,FCH3COOHoc,
                FNO3oc,FNO2oc,FN2oc,FH2SO4oc,FH2Soc,FX0Doc,time1,T0,T,fact,pCO2,Volc,tmax,
                t_exclu_m,FTnames,FTlist,t_intro_list,Delta,alphaC,Ccycle)

        time_vec.append(time2)

        ################################################
        # OSCILLATION ATTENUATIONS:
        ############

        [Hatm0,Catm0,Gatm0,COatm0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0,T] = \
        oscillation_PB(Hatm0,Catm0,Gatm0,COatm0,T_0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,
                       X0Doc0,
                       Hatm1,Catm1,Gatm1,COatm1,T1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,
                       X0Doc1,
                       Hatm2,Catm2,Gatm2,COatm2,T2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,
                       X0Doc2,
                       time0,time1,time2,tmax)

        alphaC = exp(9345.17/T-167.8108+23.3585*log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35.0)

        count += 1
        
        if count >= 10 or time2 >= tmax:
            for i in FTnames:
                if FTlist[i] is True:
                    NCT                = np.array(Bio[i][0])
                    X0T                = np.array(Bio[i][1])
                    Biomass_list[i][0] = Biomass_list[i][0] + [sum(NCT[index:]*X0T[index:]*diff(t[index-1:]))
                                                               /sum(diff(t[index-1:]))]
                    Biomass_list[i][1] = Biomass_list[i][1] + [time_vec[-1]]
                    init_FT_list[i]    = np.array([NCT[-1],X0T[-1]])
            time_atm.append(time2)
            HatmT.append(Hatm0)
            CatmT.append(Catm0)
            GatmT.append(Gatm0)
            COatmT.append(COatm0)
            HocT.append(Hoc0)
            CocT.append(Coc0)
            NocT.append(Noc0)
            GocT.append(Goc0)
            COocT.append(COoc0)
            CH3COOHocT.append(CH3COOHoc0)
            NO3ocT.append(NO3oc0)
            NO2ocT.append(NO2oc0)
            N2ocT.append(N2oc0)
            H2SO4ocT.append(H2SO4oc0)
            H2SocT.append(H2Soc0)
            X0DocT.append(X0Doc0)
            CoT.append(Co)
            CpT.append(Cp)
            ALKoT.append(ALKo)
            ALKpT.append(ALKp)
            concCaoT.append(concCao)
            concCapT.append(concCap)
            T_vec.append(T)

            CH3COOH_vec.append(sum(CH3COOHT[index:] * diff(t[index-1:]))/sum(diff(t[index-1:])))
            C_vec.append(sum(CT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:])))
            Co_vec.append(Co)
            NO3_vec.append(sum(NO3T[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))
            NO2_vec.append(sum(NO2T[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))
            H2SO4_vec.append(sum(H2SO4T[index:]     * diff(t[index-1:]))/sum(diff(t[index-1:])))
            H2S_vec.append(sum(H2ST[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))
            CO_vec.append(sum(COT[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:])))
            N2_vec.append(sum(N2T[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:])))
            X0D_vec.append(sum(X0DT[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))

            FHa_vec.append(FHa*S)
            FCa_vec.append(FCa*S)
            FGa_vec.append(FGa*S)
            FCOa_vec.append(FCOa*S)
            FHoc_vec.append(FHoc)
            FCoc_vec.append(FCoc)
            FCoc_vec.append(FNoc)
            FGoc_vec.append(FGoc)
            FCOoc_vec.append(FCOoc)
            FCH3COOHoc_vec.append(FCH3COOHoc)
            FNO3oc_vec.append(FNO3oc)
            FNO2oc_vec.append(FNO2oc)
            FN2oc_vec.append(FN2oc)
            FH2SO4oc_vec.append(FH2SO4oc)
            FH2Soc_vec.append(FH2Soc)
            FX0Doc_vec.append(FX0Doc)
            NPP_vec.append(NPP)
            for i in FTnames:
                if FTlist[i] is True:
                    NPP_vec_list[i].append(NPP_list[i])
            count = 0

        Slice = 1    
        if time2/tmax * 1000 > count_dump:
            with open('dump','wb') as f:
                pickle.dump([Hatm,HatmT[::Slice],CatmT[::Slice],GatmT[::Slice],COatmT[::Slice],time_atm[::Slice],HocT[::Slice],
                             CocT[::Slice],NocT[::Slice],GocT[::Slice],COocT[::Slice],CH3COOHocT[::Slice],NO3ocT[::Slice],
                             NO2ocT[::Slice],N2ocT[::Slice],H2SO4ocT[::Slice],H2SocT[::Slice],X0DocT[::Slice],CoT[::Slice],
                             CpT[::Slice],ALKoT[::Slice],ALKpT[::Slice],concCaoT[::Slice],concCapT[::Slice],FHa_vec[::Slice],
                             FCa_vec[::Slice],FGa_vec[::Slice],FHoc_vec[::Slice],FCoc_vec[::Slice],FNoc_vec[::Slice],
                             FGoc_vec[::Slice],FCOoc_vec[::Slice],FCH3COOHoc_vec[::Slice],FNO3oc_vec[::Slice],
                             FNO2oc_vec[::Slice],FN2oc_vec[::Slice],FH2SO4oc_vec[::Slice],FH2Soc_vec[::Slice],
                             FX0Doc_vec[::Slice],CH3COOH_vec[::Slice],C_vec[::Slice],Co_vec[::Slice],NO3_vec[::Slice],
                             NO2_vec[::Slice],N2_vec[::Slice],H2SO4_vec[::Slice],H2S_vec[::Slice],CO_vec[::Slice],
                             X0D_vec[::Slice],NPP_vec[::Slice],[NPP_vec_list[i][::Slice] for i in FTnames],time_vec[::Slice]
                             ,T_vec[::Slice],[Biomass_list[i][::Slice] for i in FTnames],FTlist],f)
            count_dump = np.ceil(time2/tmax * 1000)

        # Updates the local oceanic environment for the biological species
        Hinf       = Hatm0/ntot  * alphaH
        Cinf       = Catm0/ntot  * alphaC    # For the moment [CO2]atm. is assumed constant
        Ninf       = Noc0
        Ginf       = Gatm0/ntot  * alphaG
        COinf      = COatm0/ntot * alphaCO
        N2inf      = N2atm/ntot  * alphaN2
        Env        = np.array([Hinf,Cinf,Ninf,Ginf,COinf,CH3COOHinf,NO3inf,NO2inf,N2inf,H2SO4inf,H2Sinf,X0Dinf,T,
                      Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0])

        for i in FTnames:
            if FTlist[i] is False:
                clock_intro_list[i] += time2-time0

        # Introduction of a new functional types
        if (sum([clock_intro_list[i] >= t_intro_list[i] for i in FTnames]) > 0):

            flat = False
            time_atm.append(time2)
            HatmT.append(Hatm0)
            CatmT.append(Catm0)
            GatmT.append(Gatm0)
            COatmT.append(COatm0)
            HocT.append(Hoc0)
            CocT.append(Coc0)
            NocT.append(Noc0)
            GocT.append(Goc0)
            COocT.append(COoc0)
            CH3COOHocT.append(CH3COOHoc0)
            NO3ocT.append(NO3oc0)
            NO2ocT.append(NO2oc0)
            N2ocT.append(N2oc0)
            H2SO4ocT.append(H2SO4oc0)
            H2SocT.append(H2Soc0)
            X0DocT.append(X0Doc0)
            CoT.append(Co)
            CpT.append(Cp)
            ALKoT.append(ALKo)
            ALKpT.append(ALKp)
            concCaoT.append(concCao)
            concCapT.append(concCap)
            T_vec.append(T)

            CH3COOH_vec.append(sum(CH3COOHT[index:] * diff(t[index-1:]))/sum(diff(t[index-1:])))
            C_vec.append(sum(CT[index:]             * diff(t[index-1:]))/sum(diff(t[index-1:])))
            Co_vec.append(Co)
            NO3_vec.append(sum(NO3T[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))
            NO2_vec.append(sum(NO2T[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))
            H2SO4_vec.append(sum(H2SO4T[index:]     * diff(t[index-1:]))/sum(diff(t[index-1:])))
            H2S_vec.append(sum(H2ST[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))
            CO_vec.append(sum(COT[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:])))
            N2_vec.append(sum(N2T[index:]           * diff(t[index-1:]))/sum(diff(t[index-1:])))
            X0D_vec.append(sum(X0DT[index:]         * diff(t[index-1:]))/sum(diff(t[index-1:])))

            FHa_vec.append(FHa*S)
            FCa_vec.append(FCa*S)
            FGa_vec.append(FGa*S)
            FCOa_vec.append(FCOa*S)
            FHoc_vec.append(FHoc)
            FCoc_vec.append(FCoc)
            FCoc_vec.append(FNoc)
            FGoc_vec.append(FGoc)
            FCOoc_vec.append(FCOoc)
            FCH3COOHoc_vec.append(FCH3COOHoc)
            FNO3oc_vec.append(FNO3oc)
            FNO2oc_vec.append(FNO2oc)
            FN2oc_vec.append(FN2oc)
            FH2SO4oc_vec.append(FH2SO4oc)
            FH2Soc_vec.append(FH2Soc)
            FX0Doc_vec.append(FX0Doc)
            NPP_vec.append(NPP)
            for i in FTnames:
                if FTlist[i] is True:
                    NPP_vec_list[i].append(NPP_list[i])

        for i in FTnames:
            if clock_intro_list[i] >= t_intro_list[i]:
                FTlist[i]           = True
                init_FT_list[i]     = [1e-2,5e-11]
                clock_intro_list[i] = -tmax*10
                index_slope         = len(HatmT)
                t_intro_list[i]     = 2*tmax



        biotmax = 365 * 1

        t_clock    = libtime.time()
        prop_done  = (time2-time0)/tmax
        t_done     = t_clock - t0
        prog_speed = prop_done / t_done

        print("Progress: {0:.1e} %\t || Speed: {1:.1e} %.s-1 \t || Time: {2:.0f}".format((time2)/(tmax) * 100,100*prog_speed,t_clock-t_0), end="\r")

        time0   = time2

    #print('\n',time_atm[-1],'\t',Env,init_chem,init_meth,init_NO3metht,init_NO2metht,init_H2SO4metht,init_acet,init_acet2,init_acett,init_ferm)

    HatmT           = np.array(HatmT)
    CatmT           = np.array(CatmT)
    GatmT           = np.array(GatmT)
    COatmT          = np.array(COatmT)
    HocT            = np.array(HocT)
    CocT            = np.array(CocT)
    NocT            = np.array(NocT)
    GocT            = np.array(GocT)
    COocT           = np.array(COocT)
    CH3COOHocT      = np.array(CH3COOHocT)
    NO3ocT          = np.array(NO3ocT)
    NO2ocT          = np.array(NO2ocT)
    N2ocT           = np.array(N2ocT)
    H2SO4ocT        = np.array(H2SO4ocT)
    H2SocT          = np.array(H2SocT)
    X0DocT          = np.array(X0DocT)
    CoT             = np.array(CoT)
    CpT             = np.array(CpT)
    ALKoT           = np.array(ALKoT)
    ALKpT           = np.array(ALKpT)
    concCaoT        = np.array(concCaoT)
    concCapT        = np.array(concCapT)

    FHa_vec         = np.array(FHa_vec)
    FGa_vec         = np.array(FGa_vec)
    FHoc_vec        = np.array(FHoc_vec)
    FCoc_vec        = np.array(FCoc_vec)
    FNoc_vec        = np.array(FNoc_vec)
    FGoc_vec        = np.array(FGoc_vec)
    FCOoc_vec       = np.array(FCOoc_vec)
    FCH3COOHoc_vec  = np.array(FCH3COOHoc_vec)
    FNO3oc_vec      = np.array(FNO3oc_vec)
    FNO2oc_vec      = np.array(FNO2oc_vec)
    FN2oc_vec       = np.array(FN2oc_vec)
    FH2SO4oc_vec    = np.array(FH2SO4oc_vec)
    FH2Soc_vec      = np.array(FH2Soc_vec)
    FX0Doc_vec      = np.array(FX0Doc_vec)
    time_atm        = np.array(time_atm)
    CH3COOH_vec     = np.array(CH3COOH_vec)
    C_vec           = np.array(C_vec)
    Co_vec          = np.array(Co_vec)
    NO3_vec         = np.array(NO3_vec)
    NO2_vec         = np.array(NO2_vec)
    N2_vec          = np.array(N2_vec)
    H2SO4_vec       = np.array(H2SO4_vec)
    H2S_vec         = np.array(H2S_vec)
    CO_vec          = np.array(CO_vec)
    X0D_vec         = np.array(X0D_vec)
    NPP_vec         = np.array(NPP_vec)
    for i in FTnames:
        if FTlist[i] is True:
            NPP_vec_list[i] = np.array(NPP_vec_list[i])

    time_vec       = np.array(time_vec)
    T_vec          = np.array(T_vec)

    return(Hatm,HatmT,CatmT,GatmT,COatmT,time_atm,
           HocT,CocT,NocT,GocT,COocT,CH3COOHocT,NO3ocT,NO2ocT,N2ocT,H2SO4ocT,H2SocT,X0DocT,CoT,CpT,ALKoT,ALKpT,concCaoT,
           concCapT,FHa_vec,FCa_vec,FGa_vec,FHoc_vec,FCoc_vec,FNoc_vec,FGoc_vec,FCOoc_vec,
           FCH3COOHoc_vec,FNO3oc_vec,FNO2oc_vec,FN2oc_vec,FH2SO4oc_vec,FH2Soc_vec,FX0Doc_vec,
           CH3COOH_vec,C_vec,Co_vec,NO3_vec,NO2_vec,N2_vec,H2SO4_vec,H2S_vec,CO_vec,X0D_vec,NPP_vec,NPP_vec_list,
           time_vec,T_vec,Biomass_list,FTlist)
