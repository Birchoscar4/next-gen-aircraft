from scipy.optimize import minimize
import numpy as np
from wingbuckling import compute_wing_mass_custom

# Bounds for each design variable
bounds = [
    (0.02, 0.2),        # spar_cap_radius (m)
    (4, 20),            # num_stringers (int - to be handled inside objective)
    (0.001, 0.006),     # skin_thickness (m)
    (0.001, 0.01),      # spar_thickness (m)
]

def objective(x):
    spar_cap_radius = x[0]
    num_stringers = int(round(x[1]))  # Integer variable
    skin_thickness = x[2]
    spar_thickness = x[3]

    # Inject variables into the function by modifying global values or refactoring compute_wing_mass to accept them
    result = compute_wing_mass_custom(
        spar_cap_radius=spar_cap_radius,
        num_stringers=num_stringers,
        skin_thickness=skin_thickness,
        spar_thickness=spar_thickness
    )
    
    mass = result[6]
    buckled = result[7]

    if buckled:
        return 1e6 + mass  # Heavy penalty if panel buckled
    return mass

# Use differential evolution for mixed variable types
from scipy.optimize import differential_evolution
result = differential_evolution(objective, bounds, strategy='best1bin', popsize=15)

print("Optimized Design Variables:")
print("Spar Cap Radius:", result.x[0])
print("Num Stringers:", int(round(result.x[1])))
print("Skin Thickness:", result.x[2])
print("Spar Thickness:", result.x[3])
print("Minimum Mass:", result.fun)
