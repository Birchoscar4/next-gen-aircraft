import numpy as np
import math

##User inputs
#Mission parameters
g = 9.81 #In m/s
rho_sealevel = 1.225 #Sea level density of air in kg/m^3
safety_factor = 2 #Greater than the FAA regulation of 1.5 for contingency

#Aircraft parameters
v_stall = 55.40169231 #Aircraft stall speed for category C aircraft (140 knots) in m/s
load_factor = 2.5 #FAR maximum load factor during flight
v_a = v_stall*np.sqrt(load_factor) #Manoeuvre speed in m/s

c_wing = 6.69 #Aircraft average wing chord length

#Control surfaces
c_elevon = c_wing*0.3 #Elevon is 30% of chord
c_slats = c_wing*0.2 #Slats is 20% of chord
c_rudder = c_wing*0.28 #Rudder is 28% of chord
hinge_moment_coefficient_elevon = 0.19355535 #These are placed as design constraints from XFOIL
hinge_moment_coefficient_slats = 0.19355535 #These are placed as design constraints from XFOIL
hinge_moment_coefficient_rudder = 0.1184232 #These are placed as design constraints from XFOIL
lever_arm_elevon = 0.07 #In m, these are constraints from literature
lever_arm_slats = 0.052 #In m, these are constraints from literature
lever_arm_rudder = 0.116 #In m, these are constraints from literature

#Shared landing gear
retraction_angle_deg = 60
retraction_angle_rad = math.radians(retraction_angle_deg)
tan_retraction_angle_rad = np.tan(retraction_angle_rad)
retraction_velocity = 0.2 #in rad/s

#Main landing gear
main_landing_gear_lever_arm = 0.25 #in 
main_landing_gear_mass = 946.56036 #in kg (could use the variable from Rory's CSV file?)
main_landing_gear_cofg_arm = 0.837 #m

#Nose landing gear
nose_landing_gear_lever_arm = 0.25 #in 
nose_landing_gear_mass = 169.486 #in kg (could use the variable from Rory's CSV file?)
nose_landing_gear_cofg_arm = 0.84 #m

#Cargo doors
hinge_quantity = 2
door_length = 1.82 #in m
door_width = 1.19 #in m
door_thickness = 0.1 #in m
door_material_density = 2700 #in kg/m^3
door_lever_arm = 0.25 #in m

#Brakes
total_thrust = 254621.4507 #in N, thrust required to hold the aircraft at full engine run up before takeoff
tyre_radius_ratio = 2 
brake_radius_ratio = 1
mu_brake = 0.5 #braking coefficient
main_landing_gear_quantity = 2

#Quantity of EHAs required
quantity_elevon_EHA = 8 #2 inboard and 2 outboard elevons for left and right side of aircraft, duplicated for redundancy
quantity_slats_EHA = 6 #3 per side of the aircraft, surface redundancy implemented as not flight critical
quantity_rudder_EHA = 4 #2 for each rudder for redundancy
quantity_main_landing_gear_EHA = 4 #2 for each landing gear for redundancy as flight critical
quantity_nose_landing_gear_EHA = 2 #2 for redundancy as flight critical
quantity_cargo_doors_EHA = 2 #2 for redundancy
quantity_brakes_EHA = 12 #2 for each set of wheels, 2 landing gears and triple bogey configuration gives 6 sets of wheels

##Control surface calculations
#Hinge moment evaluations
q_times_c_squared_elevon = 0.5*rho_sealevel*v_a**2*c_elevon**2
q_times_c_squared_slats = 0.5*rho_sealevel*v_a**2*c_slats**2
q_times_c_squared_rudder = 0.5*rho_sealevel*v_a**2*c_rudder**2

elevon_max_hinge_moment = q_times_c_squared_elevon*hinge_moment_coefficient_elevon
slats_max_hinge_moment = q_times_c_squared_slats*hinge_moment_coefficient_slats
rudder_max_hinge_moment = q_times_c_squared_rudder*hinge_moment_coefficient_rudder

force_elevon_EHA = elevon_max_hinge_moment/lever_arm_elevon
force_slats_EHA = slats_max_hinge_moment/lever_arm_slats
force_rudder_EHA = rudder_max_hinge_moment/lever_arm_rudder

##Landing gear force required calculations
#Main landing gear 
force_main_landing_gear_EHA = (1/main_landing_gear_lever_arm)*(2*main_landing_gear_cofg_arm*main_landing_gear_mass*g*tan_retraction_angle_rad+main_landing_gear_cofg_arm**2*main_landing_gear_mass*retraction_velocity)
#Nose landing gear
force_nose_landing_gear_EHA = (1/nose_landing_gear_lever_arm)*(2*nose_landing_gear_cofg_arm*nose_landing_gear_mass*g*tan_retraction_angle_rad+nose_landing_gear_cofg_arm**2*nose_landing_gear_mass*retraction_velocity)

##Cargo doors force calculation
door_volume = door_length*door_width*door_thickness 
door_mass = door_material_density*door_volume
door_weight = door_mass*g

door_force = door_weight/hinge_quantity
door_hinge_moment = door_force*door_length
force_cargo_door_EHA = door_hinge_moment/door_lever_arm

##Brakes force calculation
brake_EHA_quantity_per_gear = quantity_brakes_EHA/main_landing_gear_quantity
brake_force = total_thrust*(tyre_radius_ratio/brake_radius_ratio)
force_brakes_EHA = (brake_force/(2*mu_brake))/brake_EHA_quantity_per_gear

#Required force of EHAs
sf_force_elevon_EHA = force_elevon_EHA*safety_factor
sf_force_slats_EHA = force_slats_EHA*safety_factor
sf_force_rudder_EHA = force_rudder_EHA*safety_factor
sf_force_main_landing_gear_EHA = force_main_landing_gear_EHA*safety_factor
sf_force_nose_landing_gear_EHA = force_nose_landing_gear_EHA*safety_factor
sf_force_cargo_doors_EHA = force_cargo_door_EHA*safety_factor
sf_force_brakes_EHA = force_brakes_EHA*safety_factor

##Moog (2018) and Navarro (1997) actuators for linear scaling relationship
mass_moog = 204.12 #in kg
mass_navarro = 19.1 #in kg
force_moog = 428141.329 #in N
force_navarro = 59160 #in N
delta_y = mass_moog - mass_navarro
delta_x = force_moog - force_navarro
gradient = delta_y/delta_x #gradient of linear relationship
y_intercept = mass_moog-gradient*force_moog #y intercept of linear relationship

#Masses of actuators
mass_elevon_EHA = sf_force_elevon_EHA*gradient+y_intercept
mass_slats_EHA = sf_force_slats_EHA*gradient+y_intercept
mass_rudder_EHA = sf_force_rudder_EHA*gradient+y_intercept
mass_main_landing_gear_EHA = sf_force_main_landing_gear_EHA*gradient+y_intercept
mass_nose_landing_gear_EHA = sf_force_nose_landing_gear_EHA*gradient+y_intercept
mass_cargo_doors_EHA = sf_force_cargo_doors_EHA*gradient+y_intercept
mass_brakes_EHA = sf_force_brakes_EHA*gradient+y_intercept

#Total masses of actuators
total_mass_elevon_EHA = mass_elevon_EHA*quantity_elevon_EHA
total_mass_slats_EHA = mass_slats_EHA*quantity_slats_EHA
total_mass_rudder_EHA = mass_rudder_EHA*quantity_rudder_EHA
total_mass_main_landing_gear_EHA = mass_main_landing_gear_EHA*quantity_main_landing_gear_EHA
total_mass_nose_landing_gear_EHA = mass_nose_landing_gear_EHA*quantity_nose_landing_gear_EHA
total_mass_cargo_doors_EHA = mass_cargo_doors_EHA*quantity_cargo_doors_EHA
total_mass_brakes_EHA = mass_brakes_EHA*quantity_brakes_EHA

total_mass_EHA = total_mass_elevon_EHA + total_mass_slats_EHA + total_mass_rudder_EHA + total_mass_main_landing_gear_EHA + total_mass_nose_landing_gear_EHA + total_mass_cargo_doors_EHA + total_mass_brakes_EHA

print("Total elevon EHA mass:", total_mass_elevon_EHA, "kg with safety factor of", safety_factor)
print("Total slats EHA mass:", total_mass_slats_EHA, "kg with safety factor of", safety_factor)
print("Total rudder EHA mass:", total_mass_rudder_EHA, "kg with safety factor of", safety_factor)
print("Total main landing gear EHA mass:", total_mass_main_landing_gear_EHA, "kg with safety factor of", safety_factor)
print("Total nose landing gear EHA mass:", total_mass_nose_landing_gear_EHA, "kg with safety factor of", safety_factor)
print("Total brakes EHA mass:", total_mass_brakes_EHA, "kg with safety factor of", safety_factor)

print("Total EHA mass:", total_mass_EHA, "kg with safety factor of", safety_factor)