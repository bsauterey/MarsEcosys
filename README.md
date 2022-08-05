[![DOI](https://zenodo.org/badge/518750122.svg)](https://zenodo.org/badge/latestdoi/518750122)
[![Twitter Follow](https://img.shields.io/twitter/follow/BorisSauterey.svg?style=flat-square&logo=twitter&label=Follow)](https://twitter.com/BorisSauterey)
[![GitHub Follow](https://img.shields.io/github/followers/bsauterey.svg?style=flat-square&logo=github&label=Follow)](https://github.com/bsauterey)

## MarsEcosys: A Global Ecosystem Model for Early Mars

This repository comprises the code to run the planetary ecosystem model of Early Mars, the results output that were analysed in the Nature Astronomy paper, and the Jupyter notebook to open the corresponding files. You will find below instructions to use those. You can also contact me (boris.sauterey@biologie.ens.fr) for questions.

## Running the model:

Running the model requires python3 and several dependencies. Most of those dependencies are easy to install from pip. Cartopy and gdal can be tricky though especially if you are on mac. To run the model, in your terminal go to the directory in which you downloaded the code and run:

```bash
python3 Launcher_dyn_coupling_rho.py $name$
```

where *name* is the name you have chosen for your run. If you did not install all the dependencies, python will let you know what dependencies are missing and you can install them via pip. When everything is setup properly, the code will 

- save the parameter values used for the simulation as a pickle file (a format used in python to save in the same file objects of different types) under the name *name*.params
- integrate the dynamics of the ecosystem. Be patient! With the default configuration, the code can require as much as a few hours to reach the final timestep
- save the results as a pickle under the name *name*.result_dc.

In order to open the pickle files for the parameter values and results of the model you can run in a python script or Notebook the following lines:

```python
f = gzip.open(*name*.results_dc,'rb') 
t,y,pT,T_surfT,T_avT,rhoT,pCH4T,pH2T,pCO2T,pN2T,DiffHT,DiffGT = pickle.load(f)
p,pH2,pCO2,pCH4,a_eps,a_to,r_surf,T_av,ztot,a_T = pickle.load(open(*name*.params,'rb'))
f.close()
```

## Loading the data (i.e., model predictions) analyzed in the Nature Astronomy paper:

The procedure exactly the same as the one described just above but automatized to compile all the ~1000 simulations produced for each of the 3 scenarios explore. To do so you just have to run the python Notebook (to do so you will have to get Jupyter or similar) Load_data_sets.ipynb

That's it!

Boris Sauterey

