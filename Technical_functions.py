import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from math import *


def Delta_system(HatmT,GatmT,COatmT,HocT,CocT,NocT,GocT,COocT,CH3COOHocT,NO3ocT,NO2ocT,N2ocT,H2SO4ocT,H2SocT,X0DocT,
                 time_atm,window=10):
    flat  = False
    Delta = 0.1
    
    delta_Hatm      = sum(1 - (np.array(HatmT[-window:])
                               /np.array(HatmT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_Gatm      = sum(1 - (np.array(GatmT[-window:])
                               /np.array(GatmT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_COatm     = sum(1 - (np.array(COatmT[-window:])
                               /np.array(COatmT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_Hoc       = sum(1 - (np.array(HocT[-window:])
                               /np.array(HocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_Coc       = sum(1 - (np.array(CocT[-window:])
                               /np.array(CocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_Noc       = sum(1 - (np.array(NocT[-window:])
                               /np.array(NocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_Goc       = sum(1 - (np.array(GocT[-window:])
                               /np.array(GocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_COoc      = sum(1 - (np.array(COocT[-window:])
                               /np.array(COocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_CH3COOHoc = sum(1 - (np.array(CH3COOHocT[-window:])
                               /np.array(CH3COOHocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_NO3oc     = sum(1 - (np.array(NO3ocT[-window:])
                               /np.array(NO3ocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_NO2oc     = sum(1 - (np.array(NO2ocT[-window:])
                               /np.array(NO2ocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_N2oc      = sum(1 - (np.array(N2ocT[-window:])
                               /np.array(N2ocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_H2SO4oc   = sum(1 - (np.array(H2SO4ocT[-window:])
                               /np.array(H2SO4ocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_H2Soc     = sum(1 - (np.array(H2SocT[-window:])
                               /np.array(H2SocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])
    delta_X0Doc     = sum(1 - (np.array(X0DocT[-window:])
                               /np.array(X0DocT[-(window+1):-1])))/(time_atm[-1]-time_atm[-window])

    delta_system    = max([abs(delta_Hatm),abs(delta_Gatm),abs(delta_COatm),
                           abs(delta_Hoc),abs(delta_Coc),abs(delta_Noc),abs(delta_Goc),abs(delta_COoc),
                           abs(delta_CH3COOHoc),abs(delta_NO3oc),abs(delta_NO2oc),abs(delta_N2oc),abs(delta_H2SO4oc),
                           abs(delta_H2Soc),abs(delta_X0Doc)])

    if time_atm[-1] != time_atm[-2]:
        if delta_system < 1e2:
            flat  = True
            Delta = 0.01
    if time_atm[-1] != time_atm[-2]:
        if delta_system < 1e-5:
            flat  = True
            Delta = 0.001
            
    return(flat,Delta)



def oscillation_bio(HT,CT,NT,GT,COT,CH3COOHT,NO3T,NO2T,N2T,H2SO4T,H2ST,
                    H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,
                    H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,
                    X0DT,X0D1,X0D2,NPPT,NPP1,NPP2,dX0DT,dX0D1,dX0D2,
                    FTnames,FTlist,NCT_list,NC1_list,NC2_list,X0T_list,X01_list,X02_list,
                    t0,t1,t2):
    
    l_0 = [HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1]]
    l_1 = [H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1]
    l_2 = [H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2]

    if (len(np.where(np.sign([a-b for a,b in zip(l_1,l_0)])*np.sign([a-b for a,b in zip(l_2,l_1)]) == -1)[0]) == 0
        or
        t1 == t0
        or
        t2 == t1):
        HT.append(H2)
        CT.append(C2)
        NT.append(N2)
        GT.append(G2)
        COT.append(CO2)
        CH3COOHT.append(CH3COOH2)
        NO3T.append(NO32)
        NO2T.append(NO22)
        N2T.append(N22)
        H2SO4T.append(H2SO42)
        H2ST.append(H2S2)
        X0DT.append(X0D2)
        NPPT.append(NPP2)
        dX0DT.append(dX0D2)
        for i in FTnames:
            if FTlist[i] is True:
                NCT_list[i].append(NC2_list[i])
                X0T_list[i].append(X02_list[i])
    else:
        if abs(H1-HT[-1])/(t1-t0)+abs(H2-H1)/(t2-t1) != 0:
            HT.append((H1*abs(H1-HT[-1])/(t1-t0)+H2*abs(H2-H1)/(t2-t1))/
                      (abs(H1-HT[-1])/(t1-t0)+abs(H2-H1)/(t2-t1)))
        else: HT.append(H2)
        if abs(C1-CT[-1])/(t1-t0)+abs(C2-C1)/(t2-t1) != 0:
            CT.append((C1*abs(C1-CT[-1])/(t1-t0)+C2*abs(C2-C1)/(t2-t1))/
                      (abs(C1-CT[-1])/(t1-t0)+abs(C2-C1)/(t2-t1)))
        else: CT.append(C2)
        if abs(N1-NT[-1])/(t1-t0)+abs(N2-N1)/(t2-t1) != 0:
            NT.append((N1*abs(N1-NT[-1])/(t1-t0)+N2*abs(N2-N1)/(t2-t1))/
                      (abs(N1-NT[-1])/(t1-t0)+abs(N2-N1)/(t2-t1)))
        else: NT.append(N2)
        if abs(G1-GT[-1])/(t1-t0)+abs(G2-G1)/(t2-t1) != 0:
            GT.append((G1*abs(G1-GT[-1])/(t1-t0)+G2*abs(G2-G1)/(t2-t1))/
                      (abs(G1-GT[-1])/(t1-t0)+abs(G2-G1)/(t2-t1)))
        else: GT.append(G2)
        if abs(CO1-COT[-1])/(t1-t0)+abs(CO2-CO1)/(t2-t1) != 0:
            COT.append((CO1*abs(CO1-COT[-1])/(t1-t0)+CO2*abs(CO2-CO1)/(t2-t1))/
                       (abs(CO1-COT[-1])/(t1-t0)+abs(CO2-CO1)/(t2-t1)))
        else: COT.append(CO2)
        if abs(CH3COOH1-CH3COOHT[-1])/(t1-t0)+abs(CH3COOH2-CH3COOH1)/(t2-t1) != 0:
            CH3COOHT.append((CH3COOH1*abs(CH3COOH1-CH3COOHT[-1])/(t1-t0)+CH3COOH2*
                             abs(CH3COOH2-CH3COOH1)/(t2-t1))
                            /(abs(CH3COOH1-CH3COOHT[-1])/(t1-t0)+abs(CH3COOH2-CH3COOH1)/(t2-t1)))
        else: CH3COOHT.append(CH3COOH2)
        if abs(NO31-NO3T[-1])/(t1-t0)+abs(NO32-NO31)/(t2-t1) != 0:
            NO3T.append((NO31*abs(NO31-NO3T[-1])/(t1-t0)+NO32*abs(NO32-NO31)/(t2-t1))
                        /(abs(NO31-NO3T[-1])/(t1-t0)+abs(NO32-NO31)/(t2-t1)))
        else: NO3T.append(NO32)
        if abs(NO21-NO2T[-1])/(t1-t0)+abs(NO22-NO21)/(t2-t1) != 0:
            NO2T.append((NO21*abs(NO21-NO2T[-1])/(t1-t0)+NO22*abs(NO22-NO21)/(t2-t1))
                        /(abs(NO21-NO2T[-1])/(t1-t0)+abs(NO22-NO21)/(t2-t1)))
        else: NO2T.append(NO22)
        if abs(N21-N2T[-1])/(t1-t0)+abs(N22-N21)/(t2-t1) != 0:
            N2T.append((N21*abs(N21-N2T[-1])/(t1-t0)+N22*abs(N22-N21)/(t2-t1))
                       /(abs(N21-N2T[-1])/(t1-t0)+abs(N22-N21)/(t2-t1)))
        else: N2T.append(N22)
        if abs(H2SO41-H2SO4T[-1])/(t1-t0)+abs(H2SO42-H2SO41)/(t2-t1) != 0:
            H2SO4T.append((H2SO41*abs(H2SO41-H2SO4T[-1])/(t1-t0)+H2SO42*abs(H2SO42-H2SO41)/(t2-t1))
                          /(abs(H2SO41-H2SO4T[-1])/(t1-t0)+abs(H2SO42-H2SO41)/(t2-t1)))
        else: H2SO4T.append(H2SO42)
        if abs(H2S1-H2ST[-1])/(t1-t0)+abs(H2S2-H2S1)/(t2-t1) != 0:
            H2ST.append((H2S1*abs(H2S1-H2ST[-1])/(t1-t0)+H2S2*abs(H2S2-H2S1)/(t2-t1))
                        /(abs(H2S1-H2ST[-1])/(t1-t0)+abs(H2S2-H2S1)/(t2-t1)))
        else: H2ST.append(H2S2)
        if abs(X0D1-X0DT[-1])/(t1-t0)+abs(X0D2-X0D1)/(t2-t1) != 0:
            X0DT.append((X0D1*abs(X0D1-X0DT[-1])/(t1-t0)+X0D2*abs(X0D2-X0D1)/(t2-t1))
                        /(abs(X0D1-X0DT[-1])/(t1-t0)+abs(X0D2-X0D1)/(t2-t1)))
        else: X0DT.append(X0D2)
        for i in FTnames:
            if (FTlist[i] is True and 
                (NC1_list[i] != NCT_list[i][-1] and NC2_list[i] != NC1_list[i]) and
                (X01_list[i] != X0T_list[i][-1] and X02_list[i] != X01_list[i])):
                NCT_list[i].append((NC1_list[i]*abs(NC1_list[i]-NCT_list[i][-1])/(t1-t0)
                                    +NC2_list[i]*abs(NC2_list[i]-NC1_list[i])/(t2-t1))
                                   /(abs(NC1_list[i]-NCT_list[i][-1])/(t1-t0)
                                     +abs(NC2_list[i]-NC1_list[i])/(t2-t1)))
                if  (X01_list[i]-X0T_list[i][-1])*(X02_list[i]-X01_list[i]) != 0:
                    X0T_list[i].append((X01_list[i]*abs(X01_list[i]-X0T_list[i][-1])/(t1-t0)
                                        +X02_list[i]*abs(X02_list[i]-X01_list[i])/(t2-t1))
                                       /(abs(X01_list[i]-X0T_list[i][-1])/(t1-t0)
                                         +abs(X02_list[i]-X01_list[i])/(t2-t1)))
                else: X0T_list[i].append(X02_list[i])
            elif FTlist[i] is True:
                NCT_list[i].append(NC2_list[i])
                X0T_list[i].append(X02_list[i])
        if NPP1-NPPT[-1] != 0 or NPP2-NPP1 != 0:
            NPPT.append((NPP1*abs(NPP1-NPPT[-1])/(t1-t0)+NPP2*abs(NPP2-NPP1)/(t2-t1))
                        /(abs(NPP1-NPPT[-1])/(t1-t0)+abs(NPP2-NPP1)/(t2-t1)))
        else: NPPT.append(NPP2)
        if dX0D1-dX0DT[-1] != 0 or dX0D2-dX0D1 != 0:
            dX0DT.append((dX0D1*abs(dX0D1-dX0DT[-1])/(t1-t0)+dX0D2*abs(dX0D2-dX0D1)/(t2-t1))
                         /(abs(dX0D1-dX0DT[-1])/(t1-t0)+abs(dX0D2-dX0D1)/(t2-t1)))
        else: dX0DT.append(dX0D2)            
    
    return([HT,CT,NT,GT,COT,CH3COOHT,NO3T,NO2T,N2T,H2SO4T,H2ST,X0DT,NPPT,dX0DT,FTnames,FTlist,NCT_list,NC1_list,
           NC2_list,X0T_list,X01_list,X02_list])



def oscillation_geo(Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,Cp,ALKo
                    ,ALKp,concCao,concCap,Hatm1,Catm1,Gatm1,COatm1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1
                    ,H2SO4oc1,H2Soc1,X0Doc1,Co1,Cp1,ALKo1,ALKp1,concCao1,concCap1,Hatm2,Catm2,Gatm2,COatm2,Hoc2,Coc2,Noc2
                    ,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2,Co2,Cp2,ALKo2,ALKp2,concCao2,concCap2
                    ,t0,t1,t2):
        
    l_0 = [Hatm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc]
    l_1 = [Hatm1,Gatm1,COatm1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1]
    l_2 = [Hatm2,Gatm2,COatm2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2]
        
    if (len(np.where(np.sign([a-b for a,b in zip(l_1,l_0)])*np.sign([a-b for a,b in zip(l_2,l_1)]) == -1)[0]) == 0
        or
        t1 == t2
        or
        t0 == t1):
        
        [Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,
         NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,Cp,ALKo,ALKp,concCao,concCap] = \
        [Hatm2,Catm2,Gatm2,COatm2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,
         H2SO4oc2,H2Soc2,X0Doc2,Co2,Cp2,ALKo2,ALKp2,concCao2,concCap2]
        
    else:
        if abs(Hatm1-Hatm)/(t1-t0)+abs(Hatm2-Hatm1)/(t2-t1) != 0:
            Hatm        = ((Hatm1*abs(Hatm1-Hatm)/(t1-t0)+Hatm2*abs(Hatm2-Hatm1)/(t2-t1))
                       /(abs(Hatm1-Hatm)/(t1-t0)+abs(Hatm2-Hatm1)/(t2-t1)))
        else: Hatm      = Hatm2
        if abs(Catm1-Catm)/(t1-t0)+abs(Catm2-Catm1)/(t2-t1) != 0:
            Catm        = ((Catm1*abs(Catm1-Catm)/(t1-t0)+Catm2*abs(Catm2-Catm1)/(t2-t1))
                       /(abs(Catm1-Catm)/(t1-t0)+abs(Catm2-Catm1)/(t2-t1)))
        else: Catm      = Catm2
        if abs(Gatm1-Gatm)/(t1-t0)+abs(Gatm2-Gatm1)/(t2-t1) != 0:
            Gatm        = ((Gatm1*abs(Gatm1-Gatm)/(t1-t0)+Gatm2*abs(Gatm2-Gatm1)/(t2-t1))
                       /(abs(Gatm1-Gatm)/(t1-t0)+abs(Gatm2-Gatm1)/(t2-t1)))
        else: Gatm      = Gatm2
        if abs(COatm1-COatm)/(t1-t0)+abs(COatm2-COatm1)/(t2-t1) != 0:
            COatm       = ((COatm1*abs(COatm1-COatm)/(t1-t0)+COatm2*abs(COatm2-COatm1)/(t2-t1))
                       /(abs(COatm1-COatm)/(t1-t0)+abs(COatm2-COatm1)/(t2-t1)))
        else: COatm     = COatm2
        if abs(Hoc1-Hoc)/(t1-t0)+abs(Hoc2-Hoc1)/(t2-t1) != 0:
            Hoc         = ((Hoc1*abs(Hoc1-Hoc)/(t1-t0)+Hoc2*abs(Hoc2-Hoc1)/(t2-t1))
                       /(abs(Hoc1-Hoc)/(t1-t0)+abs(Hoc2-Hoc1)/(t2-t1)))
        else: Hoc       = Hoc2
        if abs(Coc1-Coc)/(t1-t0)+abs(Coc2-Coc1)/(t2-t1) != 0:
            Coc         = ((Coc1*abs(Coc1-Coc)/(t1-t0)+Coc2*abs(Coc2-Coc1)/(t2-t1))
                       /(abs(Coc1-Coc)/(t1-t0)+abs(Coc2-Coc1)/(t2-t1)))
        else: Coc       = Coc2
        if abs(Noc1-Noc)/(t1-t0)+abs(Noc2-Noc1)/(t2-t1) != 0:
            Noc         = ((Noc1*abs(Noc1-Noc)/(t1-t0)+Noc2*abs(Noc2-Noc1)/(t2-t1))
                       /(abs(Noc1-Noc)/(t1-t0)+abs(Noc2-Noc1)/(t2-t1)))
        else: Noc       = Noc2
        if abs(Goc1-Goc)/(t1-t0)+abs(Goc2-Goc1)/(t2-t1) != 0:
            Goc         = ((Goc1*abs(Goc1-Goc)/(t1-t0)+Goc2*abs(Goc2-Goc1)/(t2-t1))
                       /(abs(Goc1-Goc)/(t1-t0)+abs(Goc2-Goc1)/(t2-t1)))
        else: Goc       = Goc2
        if abs(COoc1-COoc)/(t1-t0)+abs(COoc2-COoc1)/(t2-t1) != 0:
            COoc        = ((COoc1*abs(COoc1-COoc)/(t1-t0)+COoc2*abs(COoc2-COoc1)/(t2-t1))
                       /(abs(COoc1-COoc)/(t1-t0)+abs(COoc2-COoc1)/(t2-t1)))
        else: COoc      = COoc2
        if abs(CH3COOHoc1-CH3COOHoc)/(t1-t0)+abs(CH3COOHoc2-CH3COOHoc1) != 0:
            CH3COOHoc   = ((CH3COOHoc1*abs(CH3COOHoc1-CH3COOHoc)/(t1-t0)+CH3COOHoc2*abs(CH3COOHoc2-CH3COOHoc1)/(t2-t1))
                       /(abs(CH3COOHoc1-CH3COOHoc)/(t1-t0)+abs(CH3COOHoc2-CH3COOHoc1)/(t2-t1)))
        else: CH3COOHoc = CH3COOHoc2
        if abs(NO3oc1-NO3oc)/(t1-t0)+abs(NO3oc2-NO3oc1)/(t2-t1) != 0:
            NO3oc       = ((NO3oc1*abs(NO3oc1-NO3oc)/(t1-t0)+NO3oc2*abs(NO3oc2-NO3oc1)/(t2-t1))
                       /(abs(NO3oc1-NO3oc)/(t1-t0)+abs(NO3oc2-NO3oc1)/(t2-t1)))
        else: NO3oc     = NO3oc2
        if abs(NO2oc1-NO2oc)/(t1-t0)+abs(NO2oc2-NO2oc1)/(t2-t1) != 0:
            NO2oc       = ((NO2oc1*abs(NO2oc1-NO2oc)/(t1-t0)+NO2oc2*abs(NO2oc2-NO2oc1)/(t2-t1))
                       /(abs(NO2oc1-NO2oc)/(t1-t0)+abs(NO2oc2-NO2oc1)/(t2-t1)))
        else: NO2oc     = NO2oc2
        if abs(N2oc1-N2oc)/(t1-t0)+abs(N2oc2-N2oc1)/(t2-t1) != 0:
            N2oc        = ((N2oc1*abs(N2oc1-N2oc)/(t1-t0)+N2oc2*abs(N2oc2-N2oc1)/(t2-t1))
                       /(abs(N2oc1-N2oc)/(t1-t0)+abs(N2oc2-N2oc1)/(t2-t1)))
        else: N2oc      = N2oc2
        if abs(H2SO4oc1-H2SO4oc)/(t1-t0)+abs(H2SO4oc2-H2SO4oc1)/(t2-t1) != 0:
            H2SO4oc     = ((H2SO4oc1*abs(H2SO4oc1-H2SO4oc)/(t1-t0)+H2SO4oc2*abs(H2SO4oc2-H2SO4oc1)/(t2-t1))
                       /(abs(H2SO4oc1-H2SO4oc)/(t1-t0)+abs(H2SO4oc2-H2SO4oc1)/(t2-t1)))
        else: H2SO4oc   = H2SO4oc2
        if abs(H2Soc1-H2Soc)/(t1-t0)+abs(H2Soc2-H2Soc1)/(t2-t1) != 0:
            H2Soc       = ((H2Soc1*abs(H2Soc1-H2Soc)/(t1-t0)+H2Soc2*abs(H2Soc2-H2Soc1)/(t2-t1))
                       /(abs(H2Soc1-H2Soc)/(t1-t0)+abs(H2Soc2-H2Soc1)/(t2-t1)))
        else: H2Soc     = H2Soc2
        if abs(X0Doc1-X0Doc)/(t1-t0)+abs(X0Doc2-X0Doc1)/(t2-t1) != 0:
            X0Doc      = ((X0Doc1*abs(X0Doc1-X0Doc)/(t1-t0)+X0Doc2*abs(X0Doc2-X0Doc1)/(t2-t1))
                       /(abs(X0Doc1-X0Doc)/(t1-t0)+abs(X0Doc2-X0Doc1)/(t2-t1)))
        else: X0Doc    = X0Doc2
        if abs(Co1-Co)/(t1-t0)+abs(Co2-Co1)/(t2-t1) != 0:
            Co      = ((Co1*abs(Co1-Co)/(t1-t0)+Co2*abs(Co2-Co1)/(t2-t1))
                       /(abs(Co1-Co)/(t1-t0)+abs(Co2-Co1)/(t2-t1)))
        else: Co    = Co2
        if abs(Cp1-Cp)/(t1-t0)+abs(Cp2-Cp1)/(t2-t1) != 0:
            Cp      = ((Cp1*abs(Cp1-Cp)/(t1-t0)+Cp2*abs(Cp2-Cp1)/(t2-t1))
                       /(abs(Cp1-Cp)/(t1-t0)+abs(Cp2-Cp1)/(t2-t1)))
        else: Cp    = Cp2
        if abs(ALKo1-ALKo)/(t1-t0)+abs(ALKo2-ALKo1)/(t2-t1) != 0:
            ALKo      = ((ALKo1*abs(ALKo1-ALKo)/(t1-t0)+ALKo2*abs(ALKo2-ALKo1)/(t2-t1))
                       /(abs(ALKo1-ALKo)/(t1-t0)+abs(ALKo2-ALKo1)/(t2-t1)))
        else: ALKo    = ALKo2
        if abs(ALKp1-ALKp)/(t1-t0)+abs(ALKp2-ALKp1)/(t2-t1) != 0:
            ALKp      = ((ALKp1*abs(ALKp1-ALKp)/(t1-t0)+ALKp2*abs(ALKp2-ALKp1)/(t2-t1))
                       /(abs(ALKp1-ALKp)/(t1-t0)+abs(ALKp2-ALKp1)/(t2-t1)))
        else: ALKp    = ALKp2
        if abs(concCao1-concCao)/(t1-t0)+abs(concCao2-concCao1)/(t2-t1) != 0:
            concCao      = ((concCao1*abs(concCao1-concCao)/(t1-t0)+concCao2*abs(concCao2-concCao1)/(t2-t1))
                       /(abs(concCao1-concCao)/(t1-t0)+abs(concCao2-concCao1)/(t2-t1)))
        else: concCao    = concCao2
        if abs(concCap1-concCap)/(t1-t0)+abs(concCap2-concCap1)/(t2-t1) != 0:
            concCap      = ((concCap1*abs(concCap1-concCap)/(t1-t0)+concCap2*abs(concCap2-concCap1)/(t2-t1))
                       /(abs(concCap1-concCap)/(t1-t0)+abs(concCap2-concCap1)/(t2-t1)))
        else: concCap    = concCap2
                
    return([Hatm,Catm,Gatm,COatm,Hoc,Coc,Noc,Goc,COoc,CH3COOHoc,NO3oc,NO2oc,N2oc,H2SO4oc,H2Soc,X0Doc,Co,Cp,ALKo
                    ,ALKp,concCao,concCap])



def oscillation_PB(Hatm0,Catm0,Gatm0,COatm0,T_0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0,
                   Hatm1,Catm1,Gatm1,COatm1,T1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1,
                   Hatm2,Catm2,Gatm2,COatm2,T2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2,
                   time0,time1,time2,tmax):
        
    l_0 = [Hatm0,Gatm0,COatm0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0]
    l_1 = [Hatm1,Gatm1,COatm1,Hoc1,Coc1,Noc1,Goc1,COoc1,CH3COOHoc1,NO3oc1,NO2oc1,N2oc1,H2SO4oc1,H2Soc1,X0Doc1]
    l_2 = [Hatm2,Gatm2,COatm2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2]
    
    deltat1 = time1 - time0
    deltat2 = time2 - time1
        
    if (time2 >= tmax
        or
        len(np.where(np.sign([a-b for a,b in zip(l_1,l_0)])*np.sign([a-b for a,b in zip(l_2,l_1)]) == -1)[0]) == 0
        or
        deltat1*deltat2 == 0):
        
        oscil = False

        [Hatm0,Catm0,Gatm0,COatm0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,
         NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0,T] = \
        [Hatm2,Catm2,Gatm2,COatm2,Hoc2,Coc2,Noc2,Goc2,COoc2,CH3COOHoc2,NO3oc2,
         NO2oc2,N2oc2,H2SO4oc2,H2Soc2,X0Doc2,T2]
                
    else:
        
        oscil = True

        if abs(Hatm1-Hatm0)/deltat1+abs(Hatm2-Hatm1)/deltat2 != 0:
            Hatm0    = ((Hatm1*abs(Hatm1-Hatm0)/deltat1 + Hatm2*abs(Hatm2-Hatm1)/deltat2)
                         /(abs(Hatm1-Hatm0)/deltat1+abs(Hatm2-Hatm1)/deltat2))
        else: Hatm0  = Hatm2
        if abs(Catm1-Catm0)/deltat1+abs(Catm2-Catm1)/deltat2 != 0:
            Catm0    = ((Catm1*abs(Catm1-Catm0)/deltat1 + Catm2*abs(Catm2-Catm1)/deltat2)
                         /(abs(Catm1-Catm0)/deltat1+abs(Catm2-Catm1)/deltat2))
        else: Catm0  = Catm2
        if abs(Gatm1-Gatm0)/deltat1+abs(Gatm2-Gatm1)/deltat2 !=0:
            Gatm0    = ((Gatm1*abs(Gatm1-Gatm0)/deltat1 + Gatm2*abs(Gatm2-Gatm1)/deltat2)
                         /(abs(Gatm1-Gatm0)/deltat1+abs(Gatm2-Gatm1)/deltat2))
        else: Gatm0  = Gatm2
        if abs(COatm1-COatm0)/deltat1+abs(COatm2-COatm1)/deltat2 !=0:
            COatm0   = ((COatm1*abs(COatm1-COatm0)/deltat1
                           + COatm2*abs(COatm2-COatm1)/deltat2)
                         /(abs(COatm1-COatm0)/deltat1+abs(COatm2-COatm1)/deltat2))
        else: COatm0 = COatm2
        if abs(Hoc1-Hoc0)/deltat1+abs(Hoc2-Hoc1)/deltat2 !=0:
            Hoc0     = ((Hoc1*abs(Hoc1-Hoc0)/deltat1
                           + Hoc2*abs(Hoc2-Hoc1)/deltat2)
                         /(abs(Hoc1-Hoc0)/deltat1+abs(Hoc2-Hoc1)/deltat2))
        else: Hoc0   = Hoc2
        if abs(Coc1-Coc0)/deltat1+abs(Coc2-Coc1)/deltat2 !=0:
            Coc0     = ((Coc1*abs(Coc1-Coc0)/deltat1
                           + Coc2*abs(Coc2-Coc1)/deltat2)
                         /(abs(Coc1-Coc0)/deltat1+abs(Coc2-Coc1)/deltat2))
        else: Coc0   = Coc2
        if abs(Noc1-Noc0)/deltat1+abs(Noc2-Noc1)/deltat2 !=0:
            Noc0     = ((Noc1*abs(Noc1-Noc0)/deltat1
                           + Noc2*abs(Noc2-Noc1)/deltat2)
                         /(abs(Noc1-Noc0)/deltat1+abs(Noc2-Noc1)/deltat2))
        else: Noc0   = Noc2
        if abs(Goc1-Goc0)/deltat1+abs(Goc2-Goc1)/deltat2 !=0:
            Goc0     = ((Goc1*abs(Goc1-Goc0)/deltat1
                           + Goc2*abs(Goc2-Goc1)/deltat2)
                         /(abs(Goc1-Goc0)/deltat1+abs(Goc2-Goc1)/deltat2))
        else: Goc0   = Goc2
        if abs(COoc1-COoc0)/deltat1+abs(COoc2-COoc1)/deltat2 !=0:
            COoc0    = ((COoc1*abs(COoc1-COoc0)/deltat1
                           + COoc2*abs(COoc2-COoc1)/deltat2)
                         /(abs(COoc1-COoc0)/deltat1+abs(COoc2-COoc1)/deltat2))
        else: COoc0  = COoc2
        if abs(CH3COOHoc1-CH3COOHoc0)/deltat1+abs(CH3COOHoc2-CH3COOHoc1)/deltat2 !=0:
            CH3COOHoc0 = ((CH3COOHoc1*abs(CH3COOHoc1-CH3COOHoc0)/deltat1
                           + CH3COOHoc2*abs(CH3COOHoc2-CH3COOHoc1)/deltat2)
                         /(abs(CH3COOHoc1-CH3COOHoc0)/deltat1+abs(CH3COOHoc2-CH3COOHoc1)/deltat2))
        else: CH3COOHoc0 = CH3COOHoc2
        if abs(NO3oc1-NO3oc0)/deltat1+abs(NO3oc2-NO3oc1)/deltat2 !=0:
            NO3oc0   = ((NO3oc1*abs(NO3oc1-NO3oc0)/deltat1
                           + NO3oc2*abs(NO3oc2-NO3oc1)/deltat2)
                         /(abs(NO3oc1-NO3oc0)/deltat1+abs(NO3oc2-NO3oc1)/deltat2))
        else: NO3oc0 = NO3oc2
        if abs(NO2oc1-NO2oc0)/deltat1+abs(NO2oc2-NO2oc1)/deltat2 !=0:
            NO2oc0   = ((NO2oc1*abs(NO2oc1-NO2oc0)/deltat1
                           + NO2oc2*abs(NO2oc2-NO2oc1)/deltat2)
                         /(abs(NO2oc1-NO2oc0)/deltat1+abs(NO2oc2-NO2oc1)/deltat2))
        else: NO2oc0 = NO2oc2
        if abs(N2oc1-N2oc0)/deltat1+abs(N2oc2-N2oc1)/deltat2 !=0:
            N2oc0    = ((N2oc1*abs(N2oc1-N2oc0)/deltat1
                           + N2oc2*abs(N2oc2-N2oc1)/deltat2)
                         /(abs(N2oc1-N2oc0)/deltat1+abs(N2oc2-N2oc1)/deltat2))
        else: N2oc0  = N2oc2
        if abs(H2SO4oc1-H2SO4oc0)/deltat1+abs(H2SO4oc2-H2SO4oc1)/deltat2 !=0:
            H2SO4oc0    = ((H2SO4oc1*abs(H2SO4oc1-H2SO4oc0)/deltat1
                           + H2SO4oc2*abs(H2SO4oc2-H2SO4oc1)/deltat2)
                         /(abs(H2SO4oc1-H2SO4oc0)/deltat1+abs(H2SO4oc2-H2SO4oc1)/deltat2))
        else: H2SO4oc0  = H2SO4oc2
        if abs(H2Soc1-H2Soc0)/deltat1+abs(H2Soc2-H2Soc1)/deltat2 !=0:
            H2Soc0    = ((H2Soc1*abs(H2Soc1-H2Soc0)/deltat1
                           + H2Soc2*abs(H2Soc2-H2Soc1)/deltat2)
                         /(abs(H2Soc1-H2Soc0)/deltat1+abs(H2Soc2-H2Soc1)/deltat2))
        else: H2Soc0  = H2Soc2
        if abs(X0Doc1-X0Doc0)/deltat1+abs(X0Doc2-X0Doc1)/deltat2 !=0:
            X0Doc0   = ((X0Doc1*abs(X0Doc1-X0Doc0)/deltat1
                           + X0Doc2*abs(X0Doc2-X0Doc1)/deltat2)
                         /(abs(X0Doc1-X0Doc0)/deltat1+abs(X0Doc2-X0Doc1)/deltat2))
        else: X0Doc0 = X0Doc2
        if abs(T1-T_0)/deltat1+abs(T2-T1)/deltat2 != 0:
            T = ((T1*abs(T1-T_0)/deltat1 + T2*abs(T2-T1)/deltat2)
                 /(abs(T1-T_0)/deltat1+abs(T2-T1)/deltat2))
        else:
            T = T2
                
    return([Hatm0,Catm0,Gatm0,COatm0,Hoc0,Coc0,Noc0,Goc0,COoc0,CH3COOHoc0,NO3oc0,NO2oc0,N2oc0,H2SO4oc0,H2Soc0,X0Doc0,T])



def diff(X):

    DeltaX = [x - X[i - 1] for i, x in enumerate(X)][1:]

    return(DeltaX)



def get_natural_cubic_spline_model(x, y, minval=None, maxval=None, n_knots=None, knots=None):
    """
    Get a natural cubic spline model for the data.

    For the knots, give (a) `knots` (as an array) or (b) minval, maxval and n_knots.

    If the knots are not directly specified, the resulting knots are equally
    space within the *interior* of (max, min).  That is, the endpoints are
    *not* included as knots.

    Parameters
    ----------
    x: np.array of float
        The input data
    y: np.array of float
        The outpur data
    minval: float 
        Minimum of interval containing the knots.
    maxval: float 
        Maximum of the interval containing the knots.
    n_knots: positive integer 
        The number of knots to create.
    knots: array or list of floats 
        The knots.

    Returns
    --------
    model: a model object
        The returned model will have following method:
        - predict(x):
            x is a numpy array. This will return the predicted y-values.
    """

    if knots:
        spline = NaturalCubicSpline(knots=knots)
    else:
        spline = NaturalCubicSpline(max=maxval, min=minval, n_knots=n_knots)

    p = Pipeline([
        ('nat_cubic', spline),
        ('regression', LinearRegression(fit_intercept=True))
    ])

    p.fit(x, y)

    return p


class AbstractSpline(BaseEstimator, TransformerMixin):
    """Base class for all spline basis expansions."""

    def __init__(self, max=None, min=None, n_knots=None, n_params=None, knots=None):
        if knots is None:
            if not n_knots:
                n_knots = self._compute_n_knots(n_params)
            knots = np.linspace(min, max, num=(n_knots + 2))[1:-1]
            max, min = np.max(knots), np.min(knots)
        self.knots = np.asarray(knots)

    @property
    def n_knots(self):
        return len(self.knots)

    def fit(self, *args, **kwargs):
        return self


class NaturalCubicSpline(AbstractSpline):
    """Apply a natural cubic basis expansion to an array.
    The features created with this basis expansion can be used to fit a
    piecewise cubic function under the constraint that the fitted curve is
    linear *outside* the range of the knots..  The fitted curve is continuously
    differentiable to the second order at all of the knots.
    This transformer can be created in two ways:
      - By specifying the maximum, minimum, and number of knots.
      - By specifying the cutpoints directly.  

    If the knots are not directly specified, the resulting knots are equally
    space within the *interior* of (max, min).  That is, the endpoints are
    *not* included as knots.
    Parameters
    ----------
    min: float 
        Minimum of interval containing the knots.
    max: float 
        Maximum of the interval containing the knots.
    n_knots: positive integer 
        The number of knots to create.
    knots: array or list of floats 
        The knots.
    """

    def _compute_n_knots(self, n_params):
        return n_params

    @property
    def n_params(self):
        return self.n_knots - 1

    def transform(self, X, **transform_params):
        X_spl = self._transform_array(X)
        if isinstance(X, pd.Series):
            col_names = self._make_names(X)
            X_spl = pd.DataFrame(X_spl, columns=col_names, index=X.index)
        return X_spl

    def _make_names(self, X):
        first_name = "{}_spline_linear".format(X.name)
        rest_names = ["{}_spline_{}".format(X.name, idx)
                      for idx in range(self.n_knots - 2)]
        return [first_name] + rest_names

    def _transform_array(self, X, **transform_params):
        X = X.squeeze()
        try:
            X_spl = np.zeros((X.shape[0], self.n_knots - 1))
        except IndexError: # For arrays with only one element
            X_spl = np.zeros((1, self.n_knots - 1))
        X_spl[:, 0] = X.squeeze()

        def d(knot_idx, x):
            def ppart(t): return np.maximum(0, t)

            def cube(t): return t*t*t
            numerator = (cube(ppart(x - self.knots[knot_idx]))
                         - cube(ppart(x - self.knots[self.n_knots - 1])))
            denominator = self.knots[self.n_knots - 1] - self.knots[knot_idx]
            return numerator / denominator

        for i in range(0, self.n_knots - 2):
            X_spl[:, i+1] = (d(i, X) - d(self.n_knots - 2, X)).squeeze()
        return X_spl
