import sys
import os
import numpy as np
from structural_wing_model import structural_wing_model
from Initial_Geom import build_initial_geometry

# Generate initial geometry for the aircraft
build_initial_geometry()

#Runs VLM code [Define VSPAero code] - Joseph add this in and run the code if it works with working vsppytools env.
from OpenVSPAeroCode_1 import aero_analysis
aero_analysis()

# Define convergence loop on wing deflection and twist
tip_def_array = [100]
counter = 0

tip_deflection_error = 1
while tip_deflection_error > 0.05:
    # Call VLM Function     # Creates pressure distribution file for the wing
    tip_deflection = structural_wing_model()
    tip_def_array.append(tip_deflection)
    tip_deflection_error = abs((tip_def_array[counter] - tip_def_array[counter-1])/tip_def_array[counter])
    print(tip_deflection_error)
    
    # Call update geometry function
    from structures_fsi import structures_fsi
    structures_fsi()

    counter += 1

print("FSI Converged")