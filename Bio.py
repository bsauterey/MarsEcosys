"""
This code is under MIT license. See the License.txt file.

This file contains the subroutines related to the biological model

Boris Sauterey
boris.sauterey@ens.fr
"""

import numpy as np
import time as libtime
from math import *
from Constants import *
from Technical_functions import *

def init_bio(FTnames,FTlist,t_intro_list,T):
    
    init_FT_list     = {'meth':[],'NO3metht':[],'NO2metht':[],'H2SO4metht':[],'acet':[],'acet2': [],'acett':[],'ferm':[],
                        'photoH2':[]}
    clock_intro_list = {'meth':0,'NO3metht':0,'NO2metht':0,'H2SO4metht':0,'acet':0,'acet2':0,'acett':0,'ferm':0,'photoH2':0}
    
    for i in FTnames:
        if (FTlist[i] is True or t_intro_list[i] == 0):
            traits_list     = Run_def_FT(i,T)
            FTlist[i]       = True
            init_FT_list[i] = np.array([1e-10,traits_list[9]])
            t_intro_list[i] = 1e100
            
    return(init_FT_list,FTlist,t_intro_list)




def Run_local_biology(Env,z,aT,init_chem,FTnames,init_FT_list,biotmax,FTlist,time,flat=False,rc_bifurcation=1,extinction=True):

    t_start      = libtime.time()
    t            = [0]
    HT           = [init_chem[0]]
    CT           = [init_chem[1]]
    NT           = [init_chem[2]]
    GT           = [init_chem[3]]
    COT          = [init_chem[4]]
    CH3COOHT     = [init_chem[5]]
    NO3T         = [init_chem[6]]
    NO2T         = [init_chem[7]]
    N2T          = [init_chem[8]]
    H2SO4T       = [init_chem[9]]
    H2ST         = [init_chem[10]]
    X0DT         = [init_chem[11]]
#    T            = Env[12]
    T            = Env[12] + aT*z                         # K
    dX0DT        = [0]
    Recycling    = [0]
    NPPT         = [0]
    NPPT_list    = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],'acett':[0],'ferm':[0],
                        'photoH2':[0]}
    traits_list   = {'meth':[],'NO3metht':[],'NO2metht':[],'H2SO4metht':[],'acet':[],'acet2': [],'acett':[],'ferm':[],
                        'photoH2':[]}
    NCT_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],'acett':[0],'ferm':[0],
                        'photoH2':[0]}
    X0T_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],'acett':[0],'ferm':[0],
                        'photoH2':[0]}
    Bio          = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],'acett':[0],'ferm':[0],
                        'photoH2':[0]}
    
    prompt = 0
    nspec  = 0
    for i in FTnames:
        if FTlist[i] is True : 
            nspec += 1
            traits_list[i] = Run_def_FT(i,T)
            NCT_list[i] = [init_FT_list[i][0]]
            X0T_list[i] = [init_FT_list[i][1]]
        
    #biotmax = 0.1 * 365        
    if (nspec >= 2):
        biotmax = 0.1 * 365
    if nspec != 0:
        while t[-1] < biotmax and nspec > 0:

            t_0 = libtime.time()

            ###################
            # K1:
            ##################
            # Vector storing the total effect all the functional types on chemical species concentration:
            tot_F        = [0,0,0,0,0,0,0,0,0,0,0,0]
            # total net primary productivity:
            dX0D         = 0
            NPP          = 0
            NPP_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],
                            'acett':[0],'ferm':[0],'photoH2':[0]}
            dNC_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],
                            'acett':[0],'ferm':[0],'photoH2':[0]}
            dX0_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],
                            'acett':[0],'ferm':[0],'photoH2':[0]}

            der_vec      = []
            vec          = []

            for i in FTnames:
                if i != 'photoH2':
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                        step_FT(HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1],
                                X0DT[-1],T,i,NCT_list[i][-1],X0T_list[i][-1],traits_list[i])
                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False
                            
                        F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NCT_list[i][-1],
                             (qcat*Catabolism[1]  + qana*Anabolism[1])  * NCT_list[i][-1],
                             (qcat*Catabolism[2]  + qana*Anabolism[2])  * NCT_list[i][-1],
                             (qcat*Catabolism[3]  + qana*Anabolism[3])  * NCT_list[i][-1],
                             (qcat*Catabolism[5]  + qana*Anabolism[5])  * NCT_list[i][-1],
                             (qcat*Catabolism[6]  + qana*Anabolism[6])  * NCT_list[i][-1],
                             (qcat*Catabolism[7]  + qana*Anabolism[7])  * NCT_list[i][-1],
                             (qcat*Catabolism[8]  + qana*Anabolism[8])  * NCT_list[i][-1],
                             (qcat*Catabolism[9]  + qana*Anabolism[9])  * NCT_list[i][-1],
                             (qcat*Catabolism[10] + qana*Anabolism[10]) * NCT_list[i][-1],
                             (qcat*Catabolism[11] + qana*Anabolism[11]) * NCT_list[i][-1],
                             (qcat*Catabolism[12] + qana*Anabolism[12]) * NCT_list[i][-1]]


                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(qana*Anabolism[1]*NCT_list[i][-1]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NCT_list[i][-1],X0T_list[i][-1]]
                else:
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                        step_FT(HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1],
                                X0DT[-1],T,i,NCT_list[i][-1],X0T_list[i][-1],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(Q_photo*Metabolism[0] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[1] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[2] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[3] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[5] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[6] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[7] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[8] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[9] )  * NCT_list[i][-1],
                             (Q_photo*Metabolism[10])  * NCT_list[i][-1],
                             (Q_photo*Metabolism[11])  * NCT_list[i][-1],
                             (Q_photo*Metabolism[12])  * NCT_list[i][-1]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(Q_photo*NCT_list[i][-1]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NCT_list[i][-1],X0T_list[i][-1]]

            nspec  = 0
            for i in FTnames:
                if FTlist[i] is True:
                    nspec += 1
            if nspec == 0:
                break
                
            dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
            Step_surf_chemistry(Env,HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],
                           NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1],X0DT[-1],tot_F,dX0D,z,T)

            der_vec     = np.array(der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net])
            vec         = np.array(vec     + [HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],
                                              NO2T[-1],N2T[-1],H2SO4T[-1],H2ST[-1],X0DT[-1]])          
            
            dt          = np.min(abs(vec[np.where(der_vec != 0)]/der_vec[np.where(der_vec != 0)])/100)
            if libtime.time() - t_start > 0.1:
                dt = dt * 10
            dt = min(max(1e-1, dt),10.)
            
            der1 = der_vec*dt
            vec1 = vec + der1/2
            #vec1 = vec * np.exp(der1/2/vec)

            vec1 = [max(vec1[i],1e-40) for i in np.arange(len(vec1))]
                
            NC1_list,X01_list,NPP1_list,H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1 = Read_vec(FTnames,FTlist,vec1)
            
            ###################
            # K2:
            ##################

            tot_F        = [0,0,0,0,0,0,0,0,0,0,0,0]

            dX0D         = 0
            NPP          = 0

            der_vec      = []
            vec          = []
            
            for i in FTnames:
                if i != 'photoH2':
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                        step_FT(H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1,T,
                                i,NC1_list[i],X01_list[i],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NC1_list[i],
                             (qcat*Catabolism[1]  + qana*Anabolism[1])  * NC1_list[i],
                             (qcat*Catabolism[2]  + qana*Anabolism[2])  * NC1_list[i],
                             (qcat*Catabolism[3]  + qana*Anabolism[3])  * NC1_list[i],
                             (qcat*Catabolism[5]  + qana*Anabolism[5])  * NC1_list[i],
                             (qcat*Catabolism[6]  + qana*Anabolism[6])  * NC1_list[i],
                             (qcat*Catabolism[7]  + qana*Anabolism[7])  * NC1_list[i],
                             (qcat*Catabolism[8]  + qana*Anabolism[8])  * NC1_list[i],
                             (qcat*Catabolism[9]  + qana*Anabolism[9])  * NC1_list[i],
                             (qcat*Catabolism[10] + qana*Anabolism[10]) * NC1_list[i],
                             (qcat*Catabolism[11] + qana*Anabolism[11]) * NC1_list[i],
                             (qcat*Catabolism[12] + qana*Anabolism[12]) * NC1_list[i]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(qana*Anabolism[1]*NC1_list[i]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NC1_list[i],X01_list[i]]
                else:
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                        step_FT(H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1,T,
                                i,NC1_list[i],X01_list[i],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(Q_photo*Metabolism[0] )  * NC1_list[i],
                             (Q_photo*Metabolism[1] )  * NC1_list[i],
                             (Q_photo*Metabolism[2] )  * NC1_list[i],
                             (Q_photo*Metabolism[3] )  * NC1_list[i],
                             (Q_photo*Metabolism[5] )  * NC1_list[i],
                             (Q_photo*Metabolism[6] )  * NC1_list[i],
                             (Q_photo*Metabolism[7] )  * NC1_list[i],
                             (Q_photo*Metabolism[8] )  * NC1_list[i],
                             (Q_photo*Metabolism[9] )  * NC1_list[i],
                             (Q_photo*Metabolism[10])  * NC1_list[i],
                             (Q_photo*Metabolism[11])  * NC1_list[i],
                             (Q_photo*Metabolism[12])  * NC1_list[i]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(Q_photo*NC1_list[i]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NC1_list[i],X01_list[i]]
            
            nspec  = 0
            for i in FTnames:
                if FTlist[i] is True:
                    nspec += 1
            if nspec == 0:
                break
                
            dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
            Step_surf_chemistry(Env,H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1,tot_F,dX0D,z,z_tot,T)

            der_vec     = np.array(der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net])
            vec         = np.array(vec     + [H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1])        
                        
            der2 = der_vec*dt
            vec2 = vec + der2/2
            #vec2 = vec * np.exp(der2/2/vec)
                
            vec2 = [max(vec2[i],1e-40) for i in np.arange(len(vec2))]

            NC2_list,X02_list,NPP2_list,H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2 = Read_vec(FTnames,FTlist,vec2)            
    
            ###################
            # K3:
            ##################

            tot_F        = [0,0,0,0,0,0,0,0,0,0,0,0]

            dX0D         = 0
            NPP          = 0

            der_vec      = []
            vec          = []
            
            for i in FTnames:
                if i != 'photoH2':
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                        step_FT(H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2,T,
                                i,NC2_list[i],X02_list[i],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NC2_list[i],
                             (qcat*Catabolism[1]  + qana*Anabolism[1])  * NC2_list[i],
                             (qcat*Catabolism[2]  + qana*Anabolism[2])  * NC2_list[i],
                             (qcat*Catabolism[3]  + qana*Anabolism[3])  * NC2_list[i],
                             (qcat*Catabolism[5]  + qana*Anabolism[5])  * NC2_list[i],
                             (qcat*Catabolism[6]  + qana*Anabolism[6])  * NC2_list[i],
                             (qcat*Catabolism[7]  + qana*Anabolism[7])  * NC2_list[i],
                             (qcat*Catabolism[8]  + qana*Anabolism[8])  * NC2_list[i],
                             (qcat*Catabolism[9]  + qana*Anabolism[9])  * NC2_list[i],
                             (qcat*Catabolism[10] + qana*Anabolism[10]) * NC2_list[i],
                             (qcat*Catabolism[11] + qana*Anabolism[11]) * NC2_list[i],
                             (qcat*Catabolism[12] + qana*Anabolism[12]) * NC2_list[i]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(qana*Anabolism[1]*NC2_list[i]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NC2_list[i],X02_list[i]]
                else:
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                        step_FT(H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2,T,
                                i,NC2_list[i],X02_list[i],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(Q_photo*Metabolism[0] )  * NC2_list[i],
                             (Q_photo*Metabolism[1] )  * NC2_list[i],
                             (Q_photo*Metabolism[2] )  * NC2_list[i],
                             (Q_photo*Metabolism[3] )  * NC2_list[i],
                             (Q_photo*Metabolism[5] )  * NC2_list[i],
                             (Q_photo*Metabolism[6] )  * NC2_list[i],
                             (Q_photo*Metabolism[7] )  * NC2_list[i],
                             (Q_photo*Metabolism[8] )  * NC2_list[i],
                             (Q_photo*Metabolism[9] )  * NC2_list[i],
                             (Q_photo*Metabolism[10])  * NC2_list[i],
                             (Q_photo*Metabolism[11])  * NC2_list[i],
                             (Q_photo*Metabolism[12])  * NC2_list[i]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(Q_photo*NC2_list[i]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NC2_list[i],X02_list[i]]
            
            nspec  = 0
            for i in FTnames:
                if FTlist[i] is True:
                    nspec += 1
            if nspec == 0:
                break
                
            dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
            Step_surf_chemistry(Env,H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2,tot_F,dX0D,z,z_tot,T)

            der_vec     = np.array(der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net])
            vec         = np.array(vec     + [H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2])        
                        
            der3 = der_vec*dt
            vec3 = vec + der3
            #vec3 = vec * np.exp(der3/vec)
            
            vec3 = [max(vec3[i],1e-40) for i in np.arange(len(vec3))]

            NC3_list,X03_list,NPP3_list,H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3 = Read_vec(FTnames,FTlist,vec3) 
            
            ###################
            # K4:
            ##################

            tot_F        = [0,0,0,0,0,0,0,0,0,0,0,0]

            dX0D         = 0
            NPP          = 0

            der_vec      = []
            vec          = []
            
            for i in FTnames:
                if i != 'photoH2':
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                        step_FT(H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3,T,
                                i,NC3_list[i],X03_list[i],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NC3_list[i],
                             (qcat*Catabolism[1]  + qana*Anabolism[1])  * NC3_list[i],
                             (qcat*Catabolism[2]  + qana*Anabolism[2])  * NC3_list[i],
                             (qcat*Catabolism[3]  + qana*Anabolism[3])  * NC3_list[i],
                             (qcat*Catabolism[5]  + qana*Anabolism[5])  * NC3_list[i],
                             (qcat*Catabolism[6]  + qana*Anabolism[6])  * NC3_list[i],
                             (qcat*Catabolism[7]  + qana*Anabolism[7])  * NC3_list[i],
                             (qcat*Catabolism[8]  + qana*Anabolism[8])  * NC3_list[i],
                             (qcat*Catabolism[9]  + qana*Anabolism[9])  * NC3_list[i],
                             (qcat*Catabolism[10] + qana*Anabolism[10]) * NC3_list[i],
                             (qcat*Catabolism[11] + qana*Anabolism[11]) * NC3_list[i],
                             (qcat*Catabolism[12] + qana*Anabolism[12]) * NC3_list[i]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(qana*Anabolism[1]*NC3_list[i]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NC3_list[i],X03_list[i]]
                else:
                    if FTlist[i] is True:
                        dNC_list[i],dX0_list[i],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                        step_FT(H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3,T,
                                i,NC3_list[i],X03_list[i],traits_list[i])

                        if NCT_list[i][-1] <= 1e-20 and dNC_list[i] < 0:
                            NCT_list[i][-1] = 1e-20
                            dNC = 0
                            if extinction is True: FTlist[i] = False

                        F = [(Q_photo*Metabolism[0] )  * NC3_list[i],
                             (Q_photo*Metabolism[1] )  * NC3_list[i],
                             (Q_photo*Metabolism[2] )  * NC3_list[i],
                             (Q_photo*Metabolism[3] )  * NC3_list[i],
                             (Q_photo*Metabolism[5] )  * NC3_list[i],
                             (Q_photo*Metabolism[6] )  * NC3_list[i],
                             (Q_photo*Metabolism[7] )  * NC3_list[i],
                             (Q_photo*Metabolism[8] )  * NC3_list[i],
                             (Q_photo*Metabolism[9] )  * NC3_list[i],
                             (Q_photo*Metabolism[10])  * NC3_list[i],
                             (Q_photo*Metabolism[11])  * NC3_list[i],
                             (Q_photo*Metabolism[12])  * NC3_list[i]]

                        tot_F       = [x+y for x,y in zip(tot_F,F)]
                        NPP_list[i] = (abs(Q_photo*NC3_list[i]))
                        NPP        += NPP_list[i]
                        dX0D       += dX0Dt
                        der_vec     = der_vec + [dNC_list[i],dX0_list[i]]
                        vec         = vec + [NC3_list[i],X03_list[i]]
            
            nspec  = 0
            for i in FTnames:
                if FTlist[i] is True:
                    nspec += 1
            if nspec == 0:
                break
                
            dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
            Step_surf_chemistry(Env,H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3,tot_F,dX0D,z,z_tot,T)

            der_vec     = np.array(der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net])
            vec         = np.array(vec     + [H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3])        
                        
            der4 = der_vec*dt
            vec4 = vec + der4
            #vec4 = vec * np.exp(der4/vec)
            
            NC4_list,X04_list,NPP4_list,H4,C4,N4,G4,CO4,CH3COOH4,NO34,NO24,N24,H2SO44,H2S4,X0D4 = Read_vec(FTnames,FTlist,vec4) 
            
            der_tot = (der1 + 2*der2 + 2*der3 + der4) / 6 / dt
            
            dt          = np.min(abs(vec[np.where(der_tot != 0)]/der_tot[np.where(der_tot != 0)])/100)
            if libtime.time() - t_start > 0.1:
                dt = dt * 10
            dt = min(max(1e-1, dt),10.)

            vec_new = vec + der_tot * dt
            #vec_new = vec * np.exp((der1 + 2*der2 + 2*der3 + der4)/6/vec)
            
            vec_new = [max(vec_new[i],1e-40) for i in np.arange(len(vec_new))]
            k = 0
            for i in FTnames:
                if FTlist[i] is True:
                    NCT_list[i].append(vec_new[k])
                    X0T_list[i].append(vec_new[k+1])
                    NPPT_list[i].append(NPP_list[i])
                    k += 2
            HT        += [vec_new[k]]
            CT        += [vec_new[k+1]]
            NT        += [vec_new[k+2]]
            GT        += [vec_new[k+3]]
            COT       += [vec_new[k+4]]
            CH3COOHT  += [vec_new[k+5]]
            NO3T      += [vec_new[k+6]]
            NO2T      += [vec_new[k+7]]
            N2T       += [vec_new[k+8]]
            H2SO4T    += [vec_new[k+9]]
            H2ST      += [vec_new[k+10]]
            X0DT      += [vec_new[k+11]]
            NPPT      += [NPP]
            
            nspec  = 0
            for i in FTnames:
                if FTlist[i] is True:
                    Bio[i] = [NCT_list[i],X0T_list[i]]
                    nspec += 1
                    
            
            t += [t[-1] + dt]
            print('{0:.3f}% || dt: {1:.2} \t'.format(t[-1]*100/biotmax,dt),end='\r')
            
    HT           = np.array(HT)
    CT           = np.array(CT)
    NT           = np.array(NT)
    GT           = np.array(GT)
    COT          = np.array(COT)
    CH3COOHT     = np.array(CH3COOHT)
    NO3T         = np.array(NO3T)
    NO2T         = np.array(NO2T)
    N2T          = np.array(N2T)
    H2SO4T       = np.array(H2SO4T)
    H2ST         = np.array(H2ST)
    X0DT         = np.array(X0DT)
    NPPT         = np.array(NPPT)
    for i in FTnames:
        if FTlist[i] is True:
            NPPT_list[i] = np.array(NPPT_list[i])
    dX0DT        = np.array(dX0DT)
        

    return(HT,CT,NT,GT,COT,CH3COOHT,NO3T,NO2T,N2T,H2SO4T,H2ST,X0DT,Bio,NPPT,NPPT_list,dX0DT,t,FTlist,Recycling)



def Run_zstruct_biology(Env,z,ztot,layers,aT,init_chem,FTnames,init_FT_list,biotmax,FTlist,time,par,
                        flat=False,rc_bifurcation=1,extinction=True):

    z_vec        = np.linspace(0,ztot,layers)
    t_start      = libtime.time()
    t            = [0]
    HT           = [[init_chem[0] for i in np.arange(layers)]]
    CT           = [[init_chem[1] for i in np.arange(layers)]]
    NT           = [[init_chem[2] for i in np.arange(layers)]]
    GT           = [[init_chem[3] for i in np.arange(layers)]]
    COT          = [[init_chem[4] for i in np.arange(layers)]]
    CH3COOHT     = [[init_chem[5] for i in np.arange(layers)]]
    NO3T         = [[init_chem[6] for i in np.arange(layers)]]
    NO2T         = [[init_chem[7] for i in np.arange(layers)]]
    N2T          = [[init_chem[8] for i in np.arange(layers)]]
    H2SO4T       = [[init_chem[9] for i in np.arange(layers)]]
    H2ST         = [[init_chem[10] for i in np.arange(layers)]]
    X0DT         = [[init_chem[11] for i in np.arange(layers)]]
#    T            = Env[12]
    T_vec        = Env[12] + aT*z_vec                         # K
    dX0DT        = [0 for i in np.arange(layers)]
    Recycling    = [0 for i in np.arange(layers)]
    NPPT         = [[0 for i in np.arange(layers)]]
    NPPT_list    = {'meth':[[0 for i in np.arange(layers)]],'NO3metht':[[0 for i in np.arange(layers)]],
                    'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                    'acet':[[0 for i in np.arange(layers)]],'acet2':[[0 for i in np.arange(layers)]],
                    'acett':[[0 for i in np.arange(layers)]],'ferm':[[0 for i in np.arange(layers)]],
                    'photoH2':[[0 for i in np.arange(layers)]]}
    traits_list  = {'meth':[[0 for i in np.arange(layers)]],'NO3metht':[[0 for i in np.arange(layers)]],
                    'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                    'acet':[[0 for i in np.arange(layers)]],'acet2':[[0 for i in np.arange(layers)]],
                    'acett':[[0 for i in np.arange(layers)]],'ferm':[[0 for i in np.arange(layers)]],
                    'photoH2':[[0 for i in np.arange(layers)]]}
    NCT_list     = {'meth':[[0 for i in np.arange(layers)]],'NO3metht':[[0 for i in np.arange(layers)]],
                    'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                    'acet':[[0 for i in np.arange(layers)]],'acet2':[[0 for i in np.arange(layers)]],
                    'acett':[[0 for i in np.arange(layers)]],'ferm':[[0 for i in np.arange(layers)]],
                    'photoH2':[[0 for i in np.arange(layers)]]}
    X0T_list     = {'meth':[[0 for i in np.arange(layers)]],'NO3metht':[[0 for i in np.arange(layers)]],
                    'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                    'acet':[[0 for i in np.arange(layers)]],'acet2':[[0 for i in np.arange(layers)]],
                    'acett':[[0 for i in np.arange(layers)]],'ferm':[[0 for i in np.arange(layers)]],
                    'photoH2':[[0 for i in np.arange(layers)]]}
    Bio          = {'meth':[[0 for i in np.arange(layers)]],'NO3metht':[[0 for i in np.arange(layers)]],
                    'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                    'acet':[[0 for i in np.arange(layers)]],'acet2':[[0 for i in np.arange(layers)]],
                    'acett':[[0 for i in np.arange(layers)]],'ferm':[[0 for i in np.arange(layers)]],
                    'photoH2':[[0 for i in np.arange(layers)]]}
    
    prompt = 0
    nspec  = 0
    for i in FTnames:
        if FTlist[i] is True : 
            nspec += 1
            traits_list[i] = [Run_def_FT(i,T_vec[k]) for k in np.arange(layers)]
            if (np.shape(init_FT_list['meth'][0]) != ()):
                NCT_list[i] = init_FT_list[i][0] 
                X0T_list[i] = init_FT_list[i][1] 
            else:
                NCT_list[i] = [[init_FT_list[i][0] for k in np.arange(layers)]]
                X0T_list[i] = [[init_FT_list[i][1]  for k in np.arange(layers)]]
                
    #biotmax = 0.1 * 365        
    if (nspec >= 2):
        biotmax = 0.1 * 365
    if nspec != 0:
        while t[-1] < biotmax and nspec > 0:

            t_0 = libtime.time()

            ###################
            # K1:
            ##################
            # Vector storing the total effect all the functional types on chemical species concentration:
            tot_F        = [[0,0,0,0,0,0,0,0,0,0,0,0] for i in np.arange(layers)]
            # total net primary productivity:
            dX0D         = 0
            NPP_list     = {'meth':[0 for i in np.arange(layers)],'NO3metht':[0 for i in np.arange(layers)],
                            'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                            'acet':[0 for i in np.arange(layers)],'acet2':[0 for i in np.arange(layers)],
                            'acett':[0 for i in np.arange(layers)],'ferm':[0 for i in np.arange(layers)],
                            'photoH2':[0 for i in np.arange(layers)]}
            dNC_list     = {'meth':[0 for i in np.arange(layers)],'NO3metht':[0 for i in np.arange(layers)],
                            'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                            'acet':[0 for i in np.arange(layers)],'acet2':[0 for i in np.arange(layers)],
                            'acett':[0 for i in np.arange(layers)],'ferm':[0 for i in np.arange(layers)],
                            'photoH2':[0 for i in np.arange(layers)]}
            dX0_list     = {'meth':[0 for i in np.arange(layers)],'NO3metht':[0 for i in np.arange(layers)],
                            'NO2metht':[[0 for i in np.arange(layers)]],'H2SO4metht':[[0 for i in np.arange(layers)]],
                            'acet':[0 for i in np.arange(layers)],'acet2':[0 for i in np.arange(layers)],
                            'acett':[0 for i in np.arange(layers)],'ferm':[0 for i in np.arange(layers)],
                            'photoH2':[0 for i in np.arange(layers)]}

            der_vec      = []
            vec          = []
            NPP          = [0 for i in np.arange(layers)]

            for z_ind in np.arange(layers):
                
                z        = z_vec[z_ind]
                T        = T_vec[z_ind]                         # K

                for i in FTnames:
                    if i != 'photoH2':
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                            step_FT(HT[-1][z_ind],CT[-1][z_ind],NT[-1][z_ind],GT[-1][z_ind],COT[-1][z_ind],
                                    CH3COOHT[-1][z_ind],NO3T[-1][z_ind],NO2T[-1][z_ind],N2T[-1][z_ind],H2SO4T[-1][z_ind],
                                    H2ST[-1][z_ind],X0DT[-1][z_ind],T,i,NCT_list[i][-1][z_ind],X0T_list[i][-1][z_ind],
                                    traits_list[i][z_ind])
                            
                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[1]  + qana*Anabolism[1])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[2]  + qana*Anabolism[2])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[3]  + qana*Anabolism[3])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[5]  + qana*Anabolism[5])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[6]  + qana*Anabolism[6])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[7]  + qana*Anabolism[7])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[8]  + qana*Anabolism[8])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[9]  + qana*Anabolism[9])  * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[10] + qana*Anabolism[10]) * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[11] + qana*Anabolism[11]) * NCT_list[i][-1][z_ind] * act,
                                 (qcat*Catabolism[12] + qana*Anabolism[12]) * NCT_list[i][-1][z_ind] * act]


                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            NPP_list[i][z_ind] = (abs(qana*Anabolism[1]*NCT_list[i][-1][z_ind]))
                            NPP[z_ind]        += NPP_list[i][z_ind]              
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NCT_list[i][-1][z_ind],X0T_list[i][-1][z_ind]]
                    else:
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                            step_FT(HT[-1][z_ind],CT[-1][z_ind],NT[-1][z_ind],GT[-1][z_ind],COT[-1][z_ind],
                                    CH3COOHT[-1][z_ind],NO3T[-1][z_ind],NO2T[-1][z_ind],N2T[-1][z_ind],H2SO4T[-1][z_ind],
                                    H2ST[-1][z_ind],X0DT[-1][z_ind],T,i,NCT_list[i][-1][z_ind],X0T_list[i][-1][z_ind],
                                    traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(Q_photo*Metabolism[0] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[1] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[2] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[3] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[5] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[6] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[7] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[8] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[9] )  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[10])  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[11])  * NCT_list[i][-1][z_ind] * act,
                                 (Q_photo*Metabolism[12])  * NCT_list[i][-1][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            NPP_list[i][z_ind] = (abs(qana*Anabolism[1]*NCT_list[i][-1][z_ind]))
                            NPP[z_ind]        += NPP_list[i][z_ind]              
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NCT_list[i][-1][z_ind],X0T_list[i][-1][z_ind]]

                dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
                Step_chemistry(Env,HT[-1],CT[-1],NT[-1],GT[-1],COT[-1],CH3COOHT[-1],NO3T[-1],NO2T[-1],N2T[-1],H2SO4T[-1],
                               H2ST[-1],X0DT[-1],tot_F[z_ind],dX0D,z_ind,T_vec,z_vec,par)

                der_vec     = der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net]
                vec         = vec     + [HT[-1][z_ind],CT[-1][z_ind],NT[-1][z_ind],GT[-1][z_ind],COT[-1][z_ind],
                                                  CH3COOHT[-1][z_ind],NO3T[-1][z_ind],NO2T[-1][z_ind],N2T[-1][z_ind],
                                                  H2SO4T[-1][z_ind],H2ST[-1][z_ind],X0DT[-1][z_ind]]          
            
            der_vec     = np.array(der_vec)
            vec         = np.array(vec)
            dt          = np.min(abs(vec[np.where(der_vec != 0)]/der_vec[np.where(der_vec != 0)])/100)
            if libtime.time() - t_start > 0.1:
                dt = dt * 10
            dt = min(max(1e-1, dt),10.)

            der1 = der_vec*dt
            vec1 = vec + der1/2
            #vec1 = vec * np.exp(der1/2/vec)

            vec1 = [max(vec1[i],1e-40) for i in np.arange(len(vec1))]

            NC1_list,X01_list,H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1 = \
            Read_vec_zstruct(FTnames,FTlist,vec1,nspec,layers)
            
            ###################
            # K2:
            ##################

            tot_F        = [[0,0,0,0,0,0,0,0,0,0,0,0] for i in np.arange(layers)]

            dX0D         = 0

            der_vec      = []
            vec          = []

            for z_ind in np.arange(layers):
                
                z        = z_vec[z_ind]
                T        = T_vec[z_ind]                         # K
            
                for i in FTnames:
                    if i != 'photoH2':
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                            step_FT(H1[z_ind],C1[z_ind],N1[z_ind],G1[z_ind],CO1[z_ind],CH3COOH1[z_ind],NO31[z_ind],NO21[z_ind],
                                    N21[z_ind],H2SO41[z_ind],H2S1[z_ind],X0D1[z_ind],T,i,NC1_list[i][z_ind],X01_list[i][z_ind],
                                    traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[1]  + qana*Anabolism[1])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[2]  + qana*Anabolism[2])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[3]  + qana*Anabolism[3])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[5]  + qana*Anabolism[5])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[6]  + qana*Anabolism[6])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[7]  + qana*Anabolism[7])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[8]  + qana*Anabolism[8])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[9]  + qana*Anabolism[9])  * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[10] + qana*Anabolism[10]) * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[11] + qana*Anabolism[11]) * NC1_list[i][z_ind] * act,
                                 (qcat*Catabolism[12] + qana*Anabolism[12]) * NC1_list[i][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NC1_list[i][z_ind],X01_list[i][z_ind]]
                    else:
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                            step_FT(H1[z_ind],C1[z_ind],N1[z_ind],G1[z_ind],CO1[z_ind],CH3COOH1[z_ind],NO31[z_ind],NO21[z_ind],
                                    N21[z_ind],H2SO41[z_ind],H2S1[z_ind_ind],X0D1[z_ind],T,i,NC1_list[i][z_ind],
                                    X01_list[i][z_ind],traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(Q_photo*Metabolism[0] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[1] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[2] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[3] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[5] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[6] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[7] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[8] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[9] )  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[10])  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[11])  * NC1_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[12])  * NC1_list[i][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NC1_list[i][z_ind],X01_list[i][z_ind]]

                nspec  = 0
                for i in FTnames:
                    if FTlist[i] is True:
                        nspec += 1
                if nspec == 0:
                    break

                dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
                Step_chemistry(Env,H1,C1,N1,G1,CO1,CH3COOH1,NO31,NO21,N21,H2SO41,H2S1,X0D1,tot_F[z_ind],dX0D,z_ind,T_vec,z_vec
                              ,par)

                der_vec     = der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net]
                vec         = vec     + [H1[z_ind],C1[z_ind],N1[z_ind],G1[z_ind],CO1[z_ind],CH3COOH1[z_ind],
                                         NO31[z_ind],NO21[z_ind],N21[z_ind],H2SO41[z_ind],H2S1[z_ind],X0D1[z_ind]]        
            
            der_vec     = np.array(der_vec)
            vec         = np.array(vec)
            der2 = der_vec*dt
            vec2 = vec + der2/2
            #vec2 = vec * np.exp(der2/2/vec)
                
            vec2 = [max(vec2[i],1e-40) for i in np.arange(len(vec2))]

            NC2_list,X02_list,H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2 = \
            Read_vec_zstruct(FTnames,FTlist,vec2,nspec,layers)            
    
            ###################
            # K3:
            ##################
    
            tot_F        = [[0,0,0,0,0,0,0,0,0,0,0,0] for i in np.arange(layers)]

            dX0D         = 0

            der_vec      = []
            vec          = []

            for z_ind in np.arange(layers):
                
                z        = z_vec[z_ind]
                T        = T_vec[z_ind]                         # K
                
                for i in FTnames:
                    if i != 'photoH2':
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                            step_FT(H2[z_ind],C2[z_ind],N2[z_ind],G2[z_ind],CO2[z_ind],CH3COOH2[z_ind],NO32[z_ind],NO22[z_ind],
                                    N22[z_ind],H2SO42[z_ind],H2S2[z_ind],X0D2[z_ind],T,i,NC2_list[i][z_ind],X02_list[i][z_ind],
                                    traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[1]  + qana*Anabolism[1])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[2]  + qana*Anabolism[2])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[3]  + qana*Anabolism[3])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[5]  + qana*Anabolism[5])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[6]  + qana*Anabolism[6])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[7]  + qana*Anabolism[7])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[8]  + qana*Anabolism[8])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[9]  + qana*Anabolism[9])  * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[10] + qana*Anabolism[10]) * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[11] + qana*Anabolism[11]) * NC2_list[i][z_ind] * act,
                                 (qcat*Catabolism[12] + qana*Anabolism[12]) * NC2_list[i][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NC2_list[i][z_ind],X02_list[i][z_ind]]
                    else:
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                            step_FT(H2[z_ind],C2[z_ind],N2[z_ind],G2[z_ind],CO2[z_ind],CH3COOH2[z_ind],NO32[z_ind],NO22[z_ind],
                                    N22[z_ind],H2SO42[z_ind],H2S2[z_ind],X0D2[z_ind],T,i,NC2_list[i][z_ind],X02_list[i][z_ind],
                                    traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0
                                if extinction is True: FTlist[i] = False

                            F = [(Q_photo*Metabolism[0] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[1] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[2] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[3] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[5] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[6] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[7] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[8] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[9] )  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[10])  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[11])  * NC2_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[12])  * NC2_list[i][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NC2_list[i][z_ind],X02_list[i][z_ind]]

                nspec  = 0
                for i in FTnames:
                    if FTlist[i] is True:
                        nspec += 1
                if nspec == 0:
                    break

                dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
                Step_chemistry(Env,H2,C2,N2,G2,CO2,CH3COOH2,NO32,NO22,N22,H2SO42,H2S2,X0D2,tot_F[z_ind],dX0D,z_ind,T_vec,z_vec
                              ,par)

                der_vec     = der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net]
                vec         = vec     + [H2[z_ind],C2[z_ind],N2[z_ind],G2[z_ind],CO2[z_ind],CH3COOH2[z_ind],NO32[z_ind],
                                         NO22[z_ind],N22[z_ind],H2SO42[z_ind],H2S2[z_ind],X0D2[z_ind]]       

            der_vec     = np.array(der_vec)
            vec         = np.array(vec)                
            der3 = der_vec*dt
            vec3 = vec + der3
            #vec3 = vec * np.exp(der3/vec)
            
            vec3 = [max(vec3[i],1e-40) for i in np.arange(len(vec3))]

            NC3_list,X03_list,H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3 = \
            Read_vec_zstruct(FTnames,FTlist,vec3,nspec,layers) 
            
            ###################
            # K4:
            ##################

            tot_F        = [[0,0,0,0,0,0,0,0,0,0,0,0] for i in np.arange(layers)]

            dX0D         = 0

            der_vec      = []
            vec          = []

            for z_ind in np.arange(layers):
                
                z        = z_vec[z_ind]
                T        = T_vec[z_ind]                         # K
            
                for i in FTnames:
                    if i != 'photoH2':
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = \
                            step_FT(H3[z_ind],C3[z_ind],N3[z_ind],G3[z_ind],CO3[z_ind],CH3COOH3[z_ind],NO33[z_ind],NO23[z_ind],
                                    N23[z_ind],H2SO43[z_ind],H2S3[z_ind],X0D3[z_ind],T,i,NC3_list[i][z_ind],X03_list[i][z_ind],
                                    traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(qcat*Catabolism[0]  + qana*Anabolism[0])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[1]  + qana*Anabolism[1])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[2]  + qana*Anabolism[2])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[3]  + qana*Anabolism[3])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[5]  + qana*Anabolism[5])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[6]  + qana*Anabolism[6])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[7]  + qana*Anabolism[7])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[8]  + qana*Anabolism[8])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[9]  + qana*Anabolism[9])  * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[10] + qana*Anabolism[10]) * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[11] + qana*Anabolism[11]) * NC3_list[i][z_ind] * act,
                                 (qcat*Catabolism[12] + qana*Anabolism[12]) * NC3_list[i][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NC3_list[i][z_ind],X03_list[i][z_ind]]
                    else:
                        if FTlist[i] is True:
                            dNC_list[i][z_ind],dX0_list[i][z_ind],dX0Dt,qmet,Q_photo,mreq,Metabolism = \
                            step_FT(H3[z_ind],C3[z_ind],N3[z_ind],G3[z_ind],CO3[z_ind],CH3COOH3[z_ind],NO33[z_ind],NO23[z_ind],
                                    N23[z_ind],H2SO43[z_ind],H2S3[z_ind],X0D3[z_ind],T,i,NC3_list[i][z_ind],X03_list[i][z_ind],
                                    traits_list[i][z_ind])

                            act = 1
                            if NCT_list[i][-1][z_ind] <= 1e-20 and dNC_list[i][z_ind] < 0:
                                NCT_list[i][-1][z_ind] = 1e-20
                                dNC = 0
                                act = 0

                            F = [(Q_photo*Metabolism[0] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[1] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[2] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[3] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[5] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[6] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[7] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[8] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[9] )  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[10])  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[11])  * NC3_list[i][z_ind] * act,
                                 (Q_photo*Metabolism[12])  * NC3_list[i][z_ind] * act]

                            tot_F[z_ind]       = [x+y for x,y in zip(tot_F[z_ind],F)]
                            dX0D              += dX0Dt
                            der_vec            = der_vec + [dNC_list[i][z_ind],dX0_list[i][z_ind]]
                            vec                = vec + [NC3_list[i][z_ind],X03_list[i][z_ind_ind]]

                nspec  = 0
                for i in FTnames:
                    if FTlist[i] is True:
                        nspec += 1
                if nspec == 0:
                    break

                dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
                Step_chemistry(Env,H3,C3,N3,G3,CO3,CH3COOH3,NO33,NO23,N23,H2SO43,H2S3,X0D3,tot_F[z_ind],dX0D,z_ind,T_vec,z_vec
                              ,par)       

                der_vec     = der_vec + [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net]
                vec         = vec     + [H3[z_ind],C3[z_ind],N3[z_ind],G3[z_ind],CO3[z_ind],CH3COOH3[z_ind],NO33[z_ind],
                                         NO23[z_ind],N23[z_ind],H2SO43[z_ind],H2S3[z_ind],X0D3[z_ind]]       

            der_vec     = np.array(der_vec)
            vec         = np.array(vec)                                
            der4 = der_vec*dt
            vec4 = vec + der4
            #vec4 = vec * np.exp(der4/vec)
            
            NC4_list,X04_list,H4,C4,N4,G4,CO4,CH3COOH4,NO34,NO24,N24,H2SO44,H2S4,X0D4 = \
            Read_vec_zstruct(FTnames,FTlist,vec4,nspec,layers) 
            
            der_tot = (der1 + 2*der2 + 2*der3 + der4) / 6 / dt
            
            dt          = np.min(abs(vec[np.where(der_tot != 0)]/der_tot[np.where(der_tot != 0)])/100)
            if libtime.time() - t_start > 0.1:
                dt = dt * 10
            dt = min(max(1e-1, dt),10.)

            vec_new = vec + der_tot * dt
            #vec_new = vec * np.exp((der1 + 2*der2 + 2*der3 + der4)/6/vec)
            
            vec_new = [max(vec_new[i],1e-40) for i in np.arange(len(vec_new))]
            
            NC_list_new,X0_list_new,H_new,C_new,N_new,G_new,CO_new,CH3COOH_new,\
            NO3_new,NO2_new,N2_new,H2SO4_new,H2S_new,X0D_new = \
            Read_vec_zstruct(FTnames,FTlist,vec_new,nspec,layers) 

            k = 0
            for i in FTnames:
                if FTlist[i] is True:
                    NCT_list[i].append(NC_list_new[i])
                    X0T_list[i].append(X0_list_new[i])
                    NPPT_list[i].append(NPP_list[i])
                    k += 2
            HT        += [H_new]
            CT        += [C_new]
            NT        += [N_new]
            GT        += [G_new]
            COT       += [CO_new]
            CH3COOHT  += [CH3COOH_new]
            NO3T      += [NO3_new]
            NO2T      += [NO2_new]
            N2T       += [N2_new]
            H2SO4T    += [H2SO4_new]
            H2ST      += [H2S_new]
            X0DT      += [X0D_new]
            NPPT      += [NPP]
            
            nspec  = 0
            for i in FTnames:
                if FTlist[i] is True:
                    Bio[i] = [NCT_list[i],X0T_list[i]]
                    nspec += 1
                    
            
            t += [t[-1] + dt]
            print('{0:.3f}% || dt: {1:.2} \t'.format(t[-1]*100/biotmax,dt),end='\r')
                            
    HT           = np.array(HT)
    CT           = np.array(CT)
    NT           = np.array(NT)
    GT           = np.array(GT)
    COT          = np.array(COT)
    CH3COOHT     = np.array(CH3COOHT)
    NO3T         = np.array(NO3T)
    NO2T         = np.array(NO2T)
    N2T          = np.array(N2T)
    H2SO4T       = np.array(H2SO4T)
    H2ST         = np.array(H2ST)
    X0DT         = np.array(X0DT)
    NPPT         = np.array(NPPT)
    for i in FTnames:
        if FTlist[i] is True:
            NPPT_list[i] = np.array(NPPT_list[i])
    dX0DT        = np.array(dX0DT)
    
    return(HT,CT,NT,GT,COT,CH3COOHT,NO3T,NO2T,N2T,H2SO4T,H2ST,X0DT,Bio,NPPT,NPPT_list,dX0DT,t,FTlist,Recycling)
    #return('bla')




def Run_def_FT(FT,T,Forced_rc=False,rc=1):

    Topt       = T                                       # Optimal temperature (K -- for the moment equal to actual
                                                         # temperature)
    if Forced_rc == False:
        if   FT == 'meth'      : rc=10**(-13.23  + 0.0431 *T)
        elif FT == 'NO3metht'  : rc=10**(-15.96  + 0.054  *T)
        elif FT == 'NO2metht'  : rc=10**(-13.053 + 0.0431 *T)
        elif FT == 'H2SO4metht': rc=10**(-13.289 + 0.0432 *T)
        elif FT == 'acet'      : rc=10**(-13.21  + 0.044  *T)
        elif FT == 'acet2'     : rc=10**(-15.11  + 0.051  *T)
        elif FT == 'acett'     : rc=10**(-12.55  + 0.042  *T)
        elif FT == 'ferm'      : rc=10**(-12.84  + 0.0431 *T)
        elif FT == 'photoH2'   : rc=10**(-13.23  + 0.0431 *T)
    else:
        rc = np.array([rc for i in np.arange(np.shape(T)[0])])

    Vc         = (4/3)*pi*rc**3                          # Cell Volume (µm3)
    Qc         = (18E-15*Vc**(0.94))                     # Cell C content (molX.Cell-1)
    qmax       = np.exp(-55.76+0.1*Topt)*Vc**0.82        # Maximum metabolic rate (d-1) (Gillooly 2001)
    mg         = 24*np.exp(-46.72+0.08*Topt)*Vc**(0.67)  # Maintenance (KJ.Cell-1) (Tijhuis 93 + Litchman 2007 -- en
                                                         # tirant sur l'exposant 0.67)
    if (np.shape(T) != ()): 
        ks         = [1e-8 for i in np.arange(len(T))]   # Metabolic half-saturation constant (mol.L-1)
        kd         = [0.5 for i in np.arange(len(T))]    # Maximum decay rate (d-1)
        mort       = [0.1 for i in np.arange(len(T))]    # Basal mortality rate (d-1)
        thresh     = Qc                                  # Inflexion point of the division function (molX.Cell-1)
        slope      = [10 for i in np.arange(len(T))]     # Slope of the division function
        gmax       = [1 for i in np.arange(len(T))]      # Maximum division rate (-d)
    else:
        ks         = 1e-8                                # Metabolic half-saturation constant (mol.L-1)
        kd         = 0.5                                 # Maximum decay rate (d-1)
        mort       = 0.1                                 # Basal mortality rate (d-1)
        thresh     = Qc                                  # Inflexion point of the division function (molX.Cell-1)
        slope      = 10                                  # Slope of the division function
        gmax       = 1                                   # Maximum division rate (-d)
    
    return ([Topt,rc,Vc,Qc,ks,qmax,mg,kd,mort,thresh,slope,gmax])



def step_FT(vec,T,FT,NC,X0,traits,p):

    #H       = max(vec[0],1e-40) * (alphaH  * p*1e-5)
    #C       = max(vec[1],1e-40) \
    #        * (np.exp(9345.17/T-167.8108+23.3585*np.log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35.0) * p*1e-5)
    #N       = max(vec[2],1e-40) 
    #G       = max(vec[3],1e-40) * (alphaG  * p*1e-5)
    #CO      = max(vec[4],1e-40) * (alphaCO * p*1e-5)
    #CH3COOH = max(vec[5],1e-40)
    #NO3     = max(vec[6],1e-40)
    #NO2     = max(vec[7],1e-40)
    #N2      = max(vec[8],1e-40) * (alphaN2  * p*1e-5)
    #H2SO4   = max(vec[9],1e-40)
    #H2S     = max(vec[10],1e-40)
    #X0D     = max(vec[11],1e-40)
    
    [H,C,G,N2] = (np.array(vec)/Av*1e3)
    
    if   FT == 'meth'      : import Methanogens as ft        ; dgcat = ft.DeltaGcat(T,H/alphaH,C/alphaC,G/alphaG)
    elif FT == 'NO3metht'  : import NO3Methanotrophs as ft   ; dgcat = ft.DeltaGcat(T,C/alphaC,G/alphaG,NO3,NO2)
    elif FT == 'NO2metht'  : import NO2Methanotrophs as ft   ; dgcat = ft.DeltaGcat(T,C/alphaC,G/alphaG,NO2,N2)
    elif FT == 'H2SO4metht': import H2SO4Methanotrophs as ft ; dgcat = ft.DeltaGcat(T,C/alphaC,G/alphaG,H2SO4,H2S)
    elif FT == 'acet'      : import Acetogens1 as ft         ; dgcat = ft.DeltaGcat(T,CO/alphaCO,C/alphaC,CH3COOH)
    elif FT == 'acet2'     : import Acetogens2 as ft         ; dgcat = ft.DeltaGcat(T,C/alphaC,H/alphaH,CH3COOH)
    elif FT == 'acett'     : import Acetotrophs as ft        ; dgcat = ft.DeltaGcat(T,C/alphaC,G/alphaG,CH3COOH)
    elif FT == 'ferm'      : import Fermentors as ft         ; dgcat = ft.DeltaGcat(T,X0D,N,H/alphaH,CH3COOH,C/alphaC)
    elif FT == 'photoH2'   : import PhotoH2 as ft
        
    Vc     = traits[2]
    Qc     = traits[3]
    ks     = traits[4]
    mg     = traits[6]
    kd     = traits[7]
    mort   = traits[8]
    thresh = traits[9]
    qmax   = traits[5]
#    qmax   = traits[5]*(10*thresh-X0)/(8*thresh)
    slope  = traits[10]
    gmax   = traits[11]

    if FT != 'photoH2':
        if dgcat < 0:
            if   FT == 'meth'      : qcat = ft.QCat(dgcat,H,C,qmax,ks)
            elif FT == 'acet'      : qcat = ft.QCat(dgcat,CO,qmax,ks)
            elif FT == 'acett'     : qcat = ft.QCat(dgcat,CH3COOH,qmax,ks)
            mreq     = ft.Mreq(mg,dgcat)                                      # Rate of cell maintenance
            decay    = ft.Decay(mreq,qcat,dgcat,kd)                           # Rate of decay
            dgana    = ft.DeltaGana(T,H/alphaH,C/alphaC,N2/alphaN2,X0/Vc)     # Energy required for the anabolic reaction
            Lam      = max(-((dgana+ft.dgdiss)/dgcat),0)                      # Metabolic coupling
            Y        = ft.Yl(Lam)                                             # Metabolic stochiometry
            slim     = ft.Slim([H,C,N2],[Y[i] for i in [0,1,9]])              # Limiting substrate
            QMet_t   = ft.QMet(dgcat,qmax,ks,slim)                            # Metabolic rate
            qana     = ft.QAna(dgcat,dgana,Lam,qcat,QMet_t,mreq,qmax,ks,slim) # Anabolic rate
            qcat     = qcat                                                   # Catabolic rates
            new_cell = ft.Gamma(thresh,slope,gmax,X0,Qc)
            dNC      = (new_cell - decay - mort)*NC
            dX0      = qana - new_cell*X0

        # if there isn't:
        else:
            qcat     = 0
            qana     = 0
            dgana    = ft.DeltaGana(T,H/alphaH,C/alphaC,N2/alphaN2,X0/Vc)     # Energy required for the anabolic reaction   
            decay    = kd
            new_cell = ft.Gamma(thresh,slope,gmax,X0,Qc)
            dNC      = (new_cell - decay - mort)*NC
            dX0      = qana - new_cell*X0
        
        dX0D         = (decay+mort)*(X0*NC)
        
        #print(' T=',T,' dgc=',dgcat,' dga=',dgana,' qc=',qcat,' qa=',qana,'\n')
        return(dNC,dX0,dX0D,qana,qcat,dgcat,ft.Catabolism,ft.Anabolism)

    else:
        slim     = ft.Slim([H,C,N2],[-2.4,-1,0.1])                       # Limiting substrate
        qmet     = ft.Monod(qmax,ks,slim)                                # Total Photosynthetic energy prod.
        mreq     = ft.Mreq(mg,T)                                         # Rate of cell maintenance
        decay    = ft.Decay(mreq,qmet,kd)                                # Rate of decay
        Q_photo  = ft.Q_photo(qmet,mreq)                                 # Biomass prod. rate
        new_cell = ft.Gamma(thresh,slope,gmax,X0,Qc)
        dNC      = (new_cell - decay - mort)*NC
        dX0      = Q_photo - new_cell*X0
        
        dX0D         = (decay+mort)*(X0*NC)
        return(dNC,dX0,dX0D,qmet,Q_photo,mreq,ft.Metabolism)



#def Step_surf_chemistry(Env,H,C,N,G,CO,CH3COOH,NO3,NO2,N2,H2SO4,H2S,X0D,tot_F,dX0D,z,z_tot,T):
def Step_surf_chemistry(Env,vec,vec_d,vec_u,tot_F,dX0D,Z,T_vec,z_vec,par):
    # Recycling of dead organic matter, assumed to be immediate until the decomposers are included

    #Hinf      = Env[0]
    #Cinf      = Env[1]
    #Ginf      = Env[3]
    #COinf     = Env[4]
    #N2inf     = Env[8]
    #Hoc       = Env[13]
    #Coc       = Env[14]
    #Noc       = Env[2]
    #Goc       = Env[16]
    #COoc      = Env[17]
    #CH3COOHoc = Env[18]
    #NO3oc     = Env[19]
    #NO2oc     = Env[20]
    #N2oc      = Env[21]
    #H2SO4oc   = Env[22]
    #H2Soc     = Env[23]
    #X0Doc     = Env[24]
#
    #r_z      = 1e-6 - 1e-6/ztot*z              # m
    #eps_z    = 0.2 * np.exp(-z/ztot)           # Ø
    #to_z     = 2.5 * np.exp(-z/(3*ztot))       # Ø
  #
    #QH       = Diff_z_X(T,r_z,eps_z,to_z,MH)
    #QC       = Diff_z_X(T,r_z,eps_z,to_z,MC)
    #QG       = Diff_z_X(T,r_z,eps_z,to_z,MG)
    #QCO      = Diff_z_X(T,r_z,eps_z,to_z,MCO)
    #QN2      = Diff_z_X(T,r_z,eps_z,to_z,MN2)
    #    
    #dH       = QH    * (Hinf-H)   + tot_F[0]
    #dC       = QC    * (Cinf-C)   + tot_F[1]
    #dG       = QG    * (Ginf-G)   + tot_F[3]
    #dCO      = QCO   * (COinf-CO) + tot_F[4]
    #dN2      = QN2   * (N2inf-N2) + tot_F[8]
    #dN       = tot_F[2]
    #dCH3COOH = Fconv * (CH3COOHoc - CH3COOH) / (MZZ*365) + tot_F[5]
    ##dCH3COOH = 0
    #dNO3     = Fconv * (NO3oc   - NO3)   / (MZZ*365) + QNO3   + tot_F[6]
    #dNO2     = Fconv * (NO2oc   - NO2)   / (MZZ*365) + QNO2   + tot_F[7]
    #dH2SO4   = Fconv * (H2SO4oc - H2SO4) / (MZZ*365) + QH2SO4 + tot_F[9]
    #dH2S     = Fconv * (H2Soc   - H2S)   / (MZZ*365) + tot_F[10]
    #dX0D_net = Fconv * (X0Doc   - X0D)   / (MZZ*365) + dX0D   + tot_F[11]
    ##dX0D_net = dX0D   + tot_F[11]
    
    H         = vec[0]
    C         = vec[1]
    G         = vec[2]
    N2        = vec[3]
    
    H_d       = vec_d[0]
    C_d       = vec_d[1]
    G_d       = vec_d[2]
    N2_d      = vec_d[3]
    
    H_u       = vec_u[0]
    C_u       = vec_u[1]
    G_u       = vec_u[2]
    N2_u      = vec_u[3]
    
    T         = T_vec[Z]
    z         = z_vec[Z]
    z_d       = z_vec[Z] + np.mean(diff(z_vec))
    T_d       = T_vec[Z] + np.mean(diff(T_vec))
    T_u       = T_vec[Z] - np.mean(diff(T_vec))

    a_eps        = par[0]
    a_to         = par[1]
    r_surf       = par[2]
    
    r_z          = r_surf - r_surf/max(z_vec)*z            # cm
    eps_z        = a_eps  * np.exp(-z/max(z_vec))          # Ø
    to_z         = a_to   * np.exp(-z/(3*max(z_vec)))      # Ø
   
    r_z_d        = r_surf - r_surf/max(z_vec)*z            # cm
    eps_z_d      = a_eps  * np.exp(-z/max(z_vec))          # Ø
    to_z_d       = a_to   * np.exp(-z/(3*max(z_vec)))      # Ø
    
    QH           = max(Diff_z_X(T,r_z,eps_z,to_z,MH),0.)  / (max(z_vec)/len(z_vec)*1e5)
    QC           = max(Diff_z_X(T,r_z,eps_z,to_z,MC),0.)  / (max(z_vec)/len(z_vec)*1e5)
    QG           = max(Diff_z_X(T,r_z,eps_z,to_z,MG),0.)  / (max(z_vec)/len(z_vec)*1e5)
    QN2          = max(Diff_z_X(T,r_z,eps_z,to_z,MN2),0.) / (max(z_vec)/len(z_vec)*1e5)
    
    QH_d         = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MH),0.)  / (max(z_vec)/len(z_vec)*1e5)
    QC_d         = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MC),0.)  / (max(z_vec)/len(z_vec)*1e5)
    QG_d         = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MG),0.)  / (max(z_vec)/len(z_vec)*1e5)
    QN2_d        = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MN2),0.) / (max(z_vec)/len(z_vec)*1e5)

    dH           = (QH  * (H_u - H) - QH_d  * (H - H_d)) / (max(z_vec)/len(z_vec)*1e5) + tot_F[0]
    dC           = (QC  * (C_u - C) - QC_d  * (C - C_d)) / (max(z_vec)/len(z_vec)*1e5) + tot_F[1]
    dG           = (QG  * (G_u - G) - QG_d  * (G - G_d)) / (max(z_vec)/len(z_vec)*1e5) + tot_F[2]
    dN2          = (QG  * (G_u - G) - QG_d  * (G - G_d)) / (max(z_vec)/len(z_vec)*1e5) + tot_F[3]

    return(dH,dC,dG,dN2)
#    return(dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net)




def Step_chemistry(Env,vec,vec_d,vec_u,tot_F,dX0D,Z,T_vec,z_vec,par,ztot,p):

    H         = vec[0]
    C         = vec[1]
    G         = vec[2]
    N2        = vec[3]

    H_d       = vec_d[0]
    C_d       = vec_d[1]
    G_d       = vec_d[2]
    N2_d      = vec_d[3]
   
    H_u       = vec_u[0]
    C_u       = vec_u[1]
    G_u       = vec_u[2]
    N2_u      = vec_u[3]
    
    z         = z_vec[Z]
    T         = T_vec[Z]
    z         = z_vec[Z]
    z_d       = z_vec[Z] + np.mean(diff(z_vec))
    T_d       = T_vec[Z] + np.mean(diff(T_vec))
    T_u       = T_vec[Z] - np.mean(diff(T_vec))
    if z == z_vec[0]:
        T_u = T
         
    a_eps        = par[0]
    a_to         = par[1]
    r_surf       = par[2]
    
    r_z          = r_surf - r_surf/ztot*z            # cm
    eps_z        = a_eps  * np.exp(-z/ztot)          # Ø
    to_z         = a_to   * np.exp(-z/(3*ztot))      # Ø
   
    r_z_d        = r_surf - r_surf/ztot*z_d            # cm
    eps_z_d      = a_eps  * np.exp(-z_d/ztot)          # Ø
    to_z_d       = a_to   * np.exp(-z_d/(3*ztot))      # Ø
    
    QH           = max(Diff_z_X(T,r_z,eps_z,to_z,MH),0.)  / (ztot/len(z_vec)*1e5)
    QC           = max(Diff_z_X(T,r_z,eps_z,to_z,MC),0.)  / (ztot/len(z_vec)*1e5)
    QG           = max(Diff_z_X(T,r_z,eps_z,to_z,MG),0.)  / (ztot/len(z_vec)*1e5)
    #QCO          = max(Diff_z_X(T,r_z,eps_z,to_z,MCO),0.) / (max(z_vec)/len(z_vec)*1e4)
    QN2          = max(Diff_z_X(T,r_z,eps_z,to_z,MN2),0.) / (ztot/len(z_vec)*1e5)
    
    QH_d         = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MH),0.)  / (ztot/len(z_vec)*1e5)
    QC_d         = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MC),0.)  / (ztot/len(z_vec)*1e5)
    QG_d         = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MG),0.)  / (ztot/len(z_vec)*1e5)
    #QCO_d        = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MCO),0.) / (max(z_vec)/len(z_vec)*1e4)
    QN2_d        = max(Diff_z_X(T_d,r_z_d,eps_z_d,to_z_d,MN2),0.) / (ztot/len(z_vec)*1e5)

    if z == max(z_vec):        
        QH_d         = 0.
        QC_d         = 0.
        QG_d         = 0.
        #QCO_d        = 0
        QN2_d        = 0.
    
    dH           = (QH  * (H_u - H)   - QH_d  * (H  - H_d))  / (ztot/len(z_vec)*1e5) + tot_F[0]
    dC           = (QC  * (C_u - C)   - QC_d  * (C  - C_d))  / (ztot/len(z_vec)*1e5) + tot_F[1]
    dG           = (QG  * (G_u - G)   - QG_d  * (G  - G_d))  / (ztot/len(z_vec)*1e5) + tot_F[2]
    dN2          = (QN2 * (N2_u - N2) - QN2_d * (N2 - N2_d)) / (ztot/len(z_vec)*1e5) + tot_F[3]

    #if z == z_vec[-4]: print(dG,'\n')
    #if z == z_vec[-3]: print(dG,QG_d,G,G_d,T_d,r_z_d,eps_z_d,to_z_d,MG,'\n')
    #print(QH  * (H_u - H),-QH_d  * (H  - H_d))
    return(dH,dC,dG,dN2)


def Diff_z_X(T_z,r_z,eps_z,to_z,M):

    #DX_z     = (eps_z*r_z)/(3*to_z) * np.sqrt((8*R*T_z)/(pi*M)) * 1e4     # cm^2/s
    DX_z     = (eps_z*r_z)/(3*to_z) * np.sqrt((8*R*1e4*T_z)/(pi*M*1e-3))   # cm^2/s
                                                                         # can also be used as influx in a dm^3/L
        
    DX_z     = DX_z * 60*60*24
        
    
    return(DX_z)




def Read_vec(FTnames,FTlist,vec):

    NC_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],
                    'acett':[0],'ferm':[0],'photoH2':[0]}
    X0_list     = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],
                    'acett':[0],'ferm':[0],'photoH2':[0]}
    NPP_list    = {'meth':[0],'NO3metht':[0],'NO2metht':[0],'H2SO4metht':[0],'acet':[0],'acet2': [0],
                    'acett':[0],'ferm':[0],'photoH2':[0]}

    k = 0
    for i in FTnames:
        if FTlist[i] is True:
            NC_list[i]   = vec[k]
            X0_list[i]   = vec[k+1]
            NPP_list[i]  = NPP_list[i]
            k += 2
    H         = vec[k]
    C         = vec[k+1]
    N         = vec[k+2]
    G         = vec[k+3]
    CO        = vec[k+4]
    CH3COOH   = vec[k+5]
    NO3       = vec[k+6]
    NO2       = vec[k+7]
    N2        = vec[k+8]
    H2SO4     = vec[k+9]
    H2S       = vec[k+10]
    X0D       = vec[k+11]
    
    return(NC_list,X0_list,NPP_list,H,C,N,G,CO,CH3COOH,NO3,NO2,N2,H2SO4,H2S,X0D)




def Read_vec_zstruct(FTnames,FTlist,vec_tot,nspec,layers):

    row_size    = int(len(vec_tot)/layers)

    H           = [0 for i in np.arange(layers)]
    C           = [0 for i in np.arange(layers)]
    N           = [0 for i in np.arange(layers)]
    G           = [0 for i in np.arange(layers)]
    CO          = [0 for i in np.arange(layers)] 
    CH3COOH     = [0 for i in np.arange(layers)]      
    NO3         = [0 for i in np.arange(layers)]     
    NO2         = [0 for i in np.arange(layers)]  
    N2          = [0 for i in np.arange(layers)] 
    H2SO4       = [0 for i in np.arange(layers)]    
    H2S         = [0 for i in np.arange(layers)]      
    X0D         = [0 for i in np.arange(layers)]  
    
    NC_list     = {'meth':[0 for i in np.arange(layers)],'NO3metht':[0 for i in np.arange(layers)],
                   'NO2metht':[0 for i in np.arange(layers)],'H2SO4metht':[0 for i in np.arange(layers)],
                   'acet':[0 for i in np.arange(layers)],'acet2':[0 for i in np.arange(layers)],
                   'acett':[0 for i in np.arange(layers)],'ferm':[0 for i in np.arange(layers)],
                   'photoH2':[0 for i in np.arange(layers)]}
    X0_list     = {'meth':[0 for i in np.arange(layers)],'NO3metht':[0 for i in np.arange(layers)],
                   'NO2metht':[0 for i in np.arange(layers)],'H2SO4metht':[0 for i in np.arange(layers)],
                   'acet':[0 for i in np.arange(layers)],'acet2':[0 for i in np.arange(layers)],
                   'acett':[0 for i in np.arange(layers)],'ferm':[0 for i in np.arange(layers)],
                   'photoH2':[0 for i in np.arange(layers)]}

    for z in np.arange(layers):
        vec = vec_tot[(row_size*z):(row_size*z)+row_size]
        k = 0
        for i in FTnames:
            if FTlist[i] is True:
                NC_list[i][z]   = vec[k]
                X0_list[i][z]   = vec[k+1]
                k += 2
        H[z]         = vec[k]
        C[z]         = vec[k+1]
        N[z]         = vec[k+2]
        G[z]         = vec[k+3]
        CO[z]        = vec[k+4]
        CH3COOH[z]   = vec[k+5]
        NO3[z]       = vec[k+6]
        NO2[z]       = vec[k+7]
        N2[z]        = vec[k+8]
        H2SO4[z]     = vec[k+9]
        H2S[z]       = vec[k+10]
        X0D[z]       = vec[k+11]
    
    return(NC_list,X0_list,H,C,N,G,CO,CH3COOH,NO3,NO2,N2,H2SO4,H2S,X0D)


