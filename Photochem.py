import numpy as np
import scipy as sc
import pandas as pd
import csv
from scipy.interpolate import LinearNDInterpolator

Grid = []

with open("Grid_no_header.csv") as csvfile:
    reader = csv.reader(csvfile,delimiter=';') # change contents to floats
    #header = next(reader)
    for row in reader: # each row is a list
        Grid.append(row)
        
Grid = np.array(Grid)
Grid = Grid.astype(float)

F_H2_lin = LinearNDInterpolator(Grid[:,(3,4)],Grid[:,-4])
F_CH4_lin = LinearNDInterpolator(Grid[:,(3,4)],Grid[:,-2])

def FH2_func(pH2,pCH4):
        
    FH2 = F_H2_lin(max(min(pCH4,0.099),1e-4),max(min(pH2,0.099),1e-4))

    return(FH2)

def FCH4_func(pH2,pCH4):
        
    FCH4 = F_CH4_lin(max(min(pCH4,0.099),1e-4),max(min(pH2,0.099),1e-4))

    return(FCH4)