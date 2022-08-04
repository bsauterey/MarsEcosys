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

    [H,C,G,N2] = (np.array(vec)/Av*1e3)
    alphaC     = np.exp(9345.17/T-167.8108+23.3585*np.log(T)+(0.023517-2.3656e-4*T+4.7036e-7*T**2)*35.0)

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
