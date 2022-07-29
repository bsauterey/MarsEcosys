import numpy as np
import scipy as sc
import pandas as pd
import csv
from scipy.interpolate import interp1d
from osgeo import gdal
from math import *


gdal_data = gdal.Open('test.tif')
gdal_band = gdal_data.GetRasterBand(1)
nodataval = gdal_band.GetNoDataValue()
data_array_Z = gdal_data.ReadAsArray().astype(np.float)
if np.any(data_array_Z == nodataval):
    data_array_Z[data_array_Z == nodataval] = np.nan
data_array_Z = np.flip(data_array_Z,0)
lon = np.linspace(-180,180,np.shape(data_array_Z)[1])
lat = np.linspace(-90,90,np.shape(data_array_Z)[0])
lat_mat = np.transpose([lat.tolist() for i in range(len(lon))])

per = []
for i in range(0,900):
    per += [(sin(lat[900+i]*pi/180)-sin(lat[900+i-1]*pi/180))/2]
per = list(reversed(per))+per

rho   = []
Tsurf = []
for T_moy in np.arange(200,320,0.5):

    data_array_T = np.empty((np.shape(data_array_Z)[0],np.shape(data_array_Z)[1]))
    for i in np.arange(np.shape(data_array_Z)[0]):
        data_array_T[i] = T_moy - pi/4*20 + 20 * cos(lat[i]*pi/180) - 2.4*data_array_Z[i]*1e-3

    weighted_per = []
    for i in range(len(lat)):
        weighted_per += [np.sum(data_array_T[i] > 273)/(len(lon))*per[i]]

    weighted_Tsurf = []
    for i in range(len(lat)):
        weighted_Tsurf += [np.mean(data_array_T[i][data_array_T[i] > 273])*weighted_per[i]/np.sum(weighted_per)]

    rho   += [np.sum(weighted_per)]
    Tsurf += [np.sum(np.array(weighted_Tsurf)[np.isnan(weighted_Tsurf) != True])]

f_rho_T_273 = interp1d(np.arange(200,320,0.5), rho, kind='linear')
f_Tsurf_T_273 = interp1d(np.arange(200,320,0.5), Tsurf, kind='linear')

rho   = []
Tsurf = []
for T_moy in np.arange(200,320,0.5):

    data_array_T = np.empty((np.shape(data_array_Z)[0],np.shape(data_array_Z)[1]))
    for i in np.arange(np.shape(data_array_Z)[0]):
        data_array_T[i] = T_moy - pi/4*20 + 20 * cos(lat[i]*pi/180) - 2.4*data_array_Z[i]*1e-3

    weighted_per = []
    for i in range(len(lat)):
        weighted_per += [np.sum(data_array_T[i] > 252)/(len(lon))*per[i]]

    weighted_Tsurf = []
    for i in range(len(lat)):
        weighted_Tsurf += [np.mean(data_array_T[i][data_array_T[i] > 252])*weighted_per[i]/np.sum(weighted_per)]

    rho   += [np.sum(weighted_per)]
    Tsurf += [np.sum(np.array(weighted_Tsurf)[np.isnan(weighted_Tsurf) != True])]


f_rho_T_252 = interp1d(np.arange(200,320,0.5), rho, kind='linear')
f_Tsurf_T_252 = interp1d(np.arange(200,320,0.5), Tsurf, kind='linear')

rho   = []
Tsurf = []
for T_moy in np.arange(200,320,0.5):

    data_array_T = np.empty((np.shape(data_array_Z)[0],np.shape(data_array_Z)[1]))
    for i in np.arange(np.shape(data_array_Z)[0]):
        data_array_T[i] = T_moy - pi/4*20 + 20 * cos(lat[i]*pi/180) - 2.4*data_array_Z[i]*1e-3

    weighted_per = []
    for i in range(len(lat)):
        weighted_per += [np.sum(data_array_T[i] > 203)/(len(lon))*per[i]]

    weighted_Tsurf = []
    for i in range(len(lat)):
        weighted_Tsurf += [np.mean(data_array_T[i][data_array_T[i] > 203])*weighted_per[i]/np.sum(weighted_per)]

    rho   += [np.sum(weighted_per)]
    Tsurf += [np.sum(np.array(weighted_Tsurf)[np.isnan(weighted_Tsurf) != True])]

f_rho_T_203   = interp1d(np.arange(200,320,0.5), rho, kind='linear')
f_Tsurf_T_203 = interp1d(np.arange(200,320,0.5), Tsurf, kind='linear')

def T_to_rho273(T):

    rho = f_rho_T_273(T)

    return(rho)

def T_to_rho252(T):

    rho = f_rho_T_252(T)

    return(rho)

def T_to_rho203(T):

    rho = f_rho_T_203(T)

    return(rho)

def T_to_Tsurf273(T):

    Tsurf = f_Tsurf_T_273(T)

    return(Tsurf)

def T_to_Tsurf252(T):

    Tsurf = f_Tsurf_T_252(T)

    return(Tsurf)

def T_to_Tsurf203(T):

    Tsurf = f_Tsurf_T_203(T)

    return(Tsurf)
