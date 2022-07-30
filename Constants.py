import numpy as np

###Environmental constants
Av         = 6e23                  # Avogadro number
R          = 8.31                  # Perfect gaz constant J.(K.mol)-1
S          = 1.4e18                # Ocean surface (cm2)
g          = 3.711                 # Gravitational constant (m.s-2)
M_CO2      = 44e-3                 # CO2 molar mass (kg.mol-1)
M_N2       = 28e-3                 # N2 molar mass (kg.mol-1)
M_CH4      = 16e-3                 # CH4 molar mass (kg.mol-1)
M_H2       = 2e-3                  # H2 molar mass (kg.mol-1)

### Physical properties of elements
# Molar mass (g.mol-1)
MH         = 1
MC         = 12
MG         = 16
MCO        = 28
MCO2       = 44
MN2        = 28
MCH3COOH   = 60

# Solubilities (mol.L-1.bar-1)
alphaH     = 7.8e-4
alphaG     = 1.4e-3
alphaCO    = 1e-3
alphaN2    = 7e-4
