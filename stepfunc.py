import numpy as np
from Bio import *
from math import isnan

def stepfunc(t,ratchet,y,layers,FTnames,FTlist,traits_list,Env,z_vec,T_vec,par,ztot,p):

    nspec  = 0
    for i in FTnames:
        if FTlist[i] is True :
            nspec += 1

    if t > ratchet:
        print('{0:.0f}'.format(t),end='\r')

    der_vec = [0 for i in np.arange(len(y))]

    y[np.where(y < 1e-100)] = 1e-100
    y[np.where(np.isnan(y) == True)] = 1e-100
    for Z in np.arange(layers):

        y_z      = y[Z::layers]
        if Z == 0:
            y_z_u    = Env[:4]
            y_z_d    = y[(Z+1)::layers]
        elif Z == layers:
            y_z_u    = y[(Z-1)::layers]
            y_z_d    = y_z
        else:
            y_z_d    = y[(Z+1)::layers]
            y_z_u    = y[(Z-1)::layers]

        z        = z_vec[Z]
        T        = T_vec[Z]
        tot_F    = [0,0,0,0]

        dNC_vec  = []
        dX0_vec  = []
        dX0D     = 0

        spec_ind = len(tot_F)-1
        for i in FTnames:
            if FTlist[i] is True:
                spec_ind += 1
                dNC,dX0,dX0Dt,qana,qcat,dgcat,Catabolism,Anabolism = step_FT(y_z[:len(tot_F)],min(T,350),i,y_z[spec_ind],
                                                                             y_z[nspec+spec_ind],
                                                                             [traits_list[i][j][Z]
                                                                              for j in np.arange(len(traits_list[i]))],p)
                dNC_vec  += [dNC]
                dX0_vec  += [dX0]
                dX0D     += dX0Dt

                F = [((qcat*Catabolism[0]  + qana*Anabolism[0])  * y_z[spec_ind]) * Av * 1e-3,
                     ((qcat*Catabolism[1]  + qana*Anabolism[1])  * y_z[spec_ind]) * Av * 1e-3,
                     #(qcat*Catabolism[2]  + qana*Anabolism[2])  * y_z[spec_ind],
                     ((qcat*Catabolism[3]  + qana*Anabolism[3])  * y_z[spec_ind]) * Av * 1e-3,
                     #(qcat*Catabolism[5]  + qana*Anabolism[5])  * y_z[spec_ind],
                     #(qcat*Catabolism[6]  + qana*Anabolism[6])  * y_z[spec_ind],
                     #(qcat*Catabolism[7]  + qana*Anabolism[7])  * y_z[spec_ind],
                     #(qcat*Catabolism[8] + qana*Anabolism[8])  * y_z[spec_ind],
                     ((qcat*Catabolism[9]  + qana*Anabolism[9])  * y_z[spec_ind]) * Av * 1e-3]
                     #(qcat*Catabolism[10] + qana*Anabolism[10]) * y_z[spec_ind],
                     #(qcat*Catabolism[11] + qana*Anabolism[11]) * y_z[spec_ind],
                     #(qcat*Catabolism[12] + qana*Anabolism[12]) * y_z[spec_ind]]

                tot_F = [x+y for x,y in zip(tot_F,F)]
                #print(t,qcat,tot_F[0]/tot_F[2]*Cons_H/Cons_G)
        #dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net = \
        dH,dC,dG,dN2 = \
        Step_chemistry(Env,y_z[:len(tot_F)],y_z_d[:len(tot_F)],y_z_u[:len(tot_F)],tot_F,dX0D,Z,T_vec,z_vec,par,ztot,p)

        #der_vec_z = [dH,dC,dN,dG,dCO,dCH3COOH,dNO3,dNO2,dN2,dH2SO4,dH2S,dX0D_net]+dNC_vec+dX0_vec
        der_vec_z = [dH,dC,dG,dN2]+dNC_vec+dX0_vec
        #der_vec_z = [dH,dC,dG,dN2]+[0]+[0]

        for i in np.arange(len(der_vec_z)):
            der_vec[(i*layers)+Z] = der_vec_z[i]

        for i in np.arange(len(y)):
            if y[i] <= 1e-99 and der_vec[i] < 0:
                der_vec[i] = 0

        for i in np.arange(400,500):
            if y[i] <= 1e-99 and der_vec[i] < 0:
                der_vec[i] = 0

    return(der_vec)









                
