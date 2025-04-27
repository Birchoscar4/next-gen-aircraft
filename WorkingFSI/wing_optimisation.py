import csv
from scipy.optimize import differential_evolution
from wing_model import compute_wing_mass

# Objective: minimize mass while avoiding buckling
def objective(x):
    spar_cap_radius, num_stringers, skin_thickness, spar_thickness = x
    num_stringers = int(round(num_stringers))
    result = compute_wing_mass(spar_cap_radius, num_stringers, skin_thickness, spar_thickness)
    mass, buckled, tip_deflection, tip_twist = result[6], result[7], result[8], result[9]
    
    penalty = 0
    if buckled == 1:
        penalty += 1e6
    if tip_deflection > D_max:
        penalty += 1e5 * (tip_deflection - D_max)  # scaled penalty
    if tip_twist > twist_max:
        penalty += 1e5 * (tip_twist - twist_max) 
    
    return mass + penalty

# Design variable bounds
bounds = [
    (0.0001, 0.08),    # spar_cap_radius
    (0, 20),         # num_stringers (integer)
    (0.0001, 0.01),  # skin_thickness
    (0.0001, 0.03)    # spar_thickness
]

for twist_max in [0.03 * i for i in range(1, 11)]:
    D_max = 0.5
    # Run the optimization
    result = differential_evolution(objective, bounds, strategy='best1bin', popsize=5, tol=0.03)
    x_opt = result.x
    x_opt[1] = int(round(x_opt[1]))  # Round stringer count

    # Re-evaluate with optimal variables to get final outputs
    final = compute_wing_mass(*x_opt)
    tip_deflection, tip_twist, mass, TE_spar = final[8], final[9], final[6], final [10]

    # Append results to CSV
    with open("optimization_results.csv", "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            tip_deflection,
            tip_twist,
            mass,
            TE_spar,
            x_opt[0],  # spar_cap_radius
            x_opt[1],  # num_stringers
            x_opt[2],  # skin_thickness
            x_opt[3]   # spar_thickness
        ])

print("Optimization complete. Results saved to optimization_results.csv")
