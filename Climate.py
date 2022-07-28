import numpy as np
import scipy as sc
import pandas as pd
import csv
from scipy.interpolate import LinearNDInterpolator

Grid = []

with open("Climate_Grid_no_header.csv") as csvfile:
    reader = csv.reader(csvfile,delimiter=';') # change contents to floats
    #header = next(reader)
    for row in reader: # each row is a list
        Grid.append(row)
        
Grid = np.array(Grid)
Grid = Grid.astype(float)

pCH4  = Grid[:,0]
pCO2  = Grid[:,1]
pH2   = Grid[:,2]
P     = Grid[:,3]
T     = Grid[:,4]

f_tot_lin = LinearNDInterpolator(Grid[:,(0,2,3)],T)

def T_func(p,pH2,pCH4):
    
    if pH2+pCH4 > 0.099:
        pH2  = pH2 /(pH2+pCH4)*0.099
        pCH4 = pCH4/(pH2+pCH4)*0.099        
    
    T = f_tot_lin(min(max(pCH4,0.0005),0.099),min(max(pH2,0.0005),0.099),min(p,2.99)) #Needs an update to integrate the (low pCH4 + low p) points 

    return(T)