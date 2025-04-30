import csv
import openvsp as vsp
from Initial_Geom import build_initial_geometry
from Initial_Geom import vlm
from cabin_mass import compute_cabin_mass
from wing_model import compute_wing_mass
from vstab_model import compute_vstab_mass
from landing_gear_sizing import landing_gear_sizing
from mass_properties import write_aircraft_data_to_csv
from structure_neg_fsi import update_geom_neg # Geometry Function from FSI
from structures_pos_fsi import update_geom_pos # Geometry Function from FSI
from vlm_fsi_neg import vlm_fsi_neg
from vlm_fsi_pos import vlm_fsi_pos
from deflected_zero_csv import neg_geom_zero
from deflected_zero_pos import pos_geom_zero
import os

build_initial_geometry()
wing_id = vsp.FindGeomsWithName("WING")  # Replace "WING" with your actual component name
wing_id = wing_id[0]

import pandas as pd

def create_zeroed_csv(output_path,rows):
    columns = ['Deflection', 'Twist']
    row_count = rows

    zero_df = pd.DataFrame(0, index=range(row_count), columns=columns)

    # Write to CSV
    zero_df.to_csv(output_path, index=False)

rows = 11
create_zeroed_csv('def_twist_pos_zero.csv',rows)
create_zeroed_csv('def_twist_neg_zero.csv',rows)



max_q_density = 0.35 #temp
max_q_velocity = 210 #temp
#Fuselage Inputs
frame_spacing = 0.2 #m
fuselage_sparcap_radius = 0.05 #m
#Wing Inputs
wing_sparcap_radius = 0.04
wing_num_stringers = 10
wing_skin_thickness = 0.002
wing_spar_thickness = 0.005
LE_spar_wing = 0.1
TE_spar_wing = 0.7
#Vstab Inputs
vstab_sparcap_radius = 0.02
vstab_num_stringers = 10
vstab_skin_thickness = 0.002
vstab_spar_thickness = 0.005
LE_spar_vstab = 0.1
TE_spar_vstab = 0.7
#Landing Gear Inputs
ground_clearance = 2
nlg_distance_from_cg = 8
mlg_distance_from_cg = 1.5
sink_speed = 3.8
z_cg = 4 # from ground

# Gust load factor parameters
C_L_a = 5.01
MTOW = 60000
rho = 1.2256
mean_chord = 5
wing_area = 275
cruise_velocity = 210
gust_velocity = 13.4

def compute_gust_loads():
    mu_g = 2*(MTOW/wing_area)/(rho*mean_chord*C_L_a)
    K_g = 0.88*mu_g/(5.3+mu_g)
    n_lim_pos = 1+(K_g*0.5*rho*cruise_velocity*gust_velocity*C_L_a)/(9.81*MTOW/wing_area)
    n_lim_neg = 1-(K_g*0.5*rho*cruise_velocity*gust_velocity*C_L_a)/(9.81*MTOW/wing_area)
    #FAR 25 Minimums
    if n_lim_pos < 2.5:
        n_lim_pos = 2.5
    if n_lim_neg > -1:
        n_lim_neg = -1
    return n_lim_pos,n_lim_neg

load_factor = [0,0]
load_factor[0], load_factor[1] = compute_gust_loads()

# Fuselage curve fit equations, returns mass properties 
cabin_mass_array = compute_cabin_mass(frame_spacing, sparcap_radius=fuselage_sparcap_radius)
#aft_body_mass_array = compute_aft_body_mass()

# Vertical Stabiliser bending analysis
vstab_mass_array = compute_vstab_mass(spar_cap_radius=vstab_sparcap_radius, num_stringers=vstab_num_stringers, skin_thickness=vstab_skin_thickness, spar_thickness=vstab_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_vstab, TE_spar=TE_spar_vstab)


### FSI LOOPS HERE
tip_def_array = []
tip_def_array.append(0)
counter = 0
# Positive load factor fsi loop
tip_deflection_error = 1
while tip_deflection_error > 0.05:
    # Call wing deflection function to produce deflection and twist distribution
    vlm_fsi_pos()
    tip_deflection = compute_wing_mass(spar_cap_radius=wing_sparcap_radius, num_stringers=wing_num_stringers, skin_thickness=wing_skin_thickness, spar_thickness=wing_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_wing, TE_spar=TE_spar_wing, file_path = "update_geom_pos_DegenGeom.slc")[8]
    tip_def_array.append(tip_deflection)
    tip_deflection_error = (tip_def_array[counter] - tip_def_array[counter-1])/tip_def_array[counter]
    file_path_delete_pos = r"C:\Users\joepe\Documents\ProjectAero\OpenVSP-3.41.1-win64-Python3.11\OpenVSP-3.41.1-win64\python\update_geom_pos_DegenGeom.slc"
    if os.path.exists(file_path_delete_pos):
        os.remove(file_path_delete_pos)
        print("File deleted successfully.")
    else:
        print("The file does not exist.")
    # Call update geometry function
    update_geom_pos(wing_id)
    vsp.WriteVSPFile("update_geom_pos.vsp3")
    print("Twist and DeltaY updated successfully.")
    print("Positive load model saved successfully as update_geom_pos.vsp3")
    # Call VLM Function     # Creates pressure distribution file for the wing
    #vlm_fsi_pos()
    counter += 1
    print (counter)
    print (tip_deflection_error)

# delete this or comment out
#vlm_fsi_pos()
#wing_mass_array = compute_wing_mass(spar_cap_radius=wing_sparcap_radius, num_stringers=wing_num_stringers, skin_thickness=wing_skin_thickness, spar_thickness=wing_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_wing, TE_spar=TE_spar_wing, file_path = "update_geom_pos_DegenGeom.slc")
#Ve_pos = wing_mass_array[12]
#alpha_e_pos = wing_mass_array[13]

neg_geom_zero(wing_id)
vsp.WriteVSPFile("update_geom_neg.vsp3")
print("Twist and DeltaY updated successfully.")
print("Negative load model saved successfully as update_geom_neg.vsp3")


tip_def_array = []
tip_def_array.append(0)
counter = 0
# Negative load factor fsi loop
tip_deflection_error = 1
while tip_deflection_error > 0.05:
    # Call wing deflection function to produce deflection and twist distribution
    vlm_fsi_neg()
    tip_deflection = compute_wing_mass(spar_cap_radius=wing_sparcap_radius, num_stringers=wing_num_stringers, skin_thickness=wing_skin_thickness, spar_thickness=wing_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_wing, TE_spar=TE_spar_wing, file_path = "update_geom_neg_DegenGeom.slc")[9]
    tip_def_array.append(tip_deflection)
    tip_deflection_error = (tip_def_array[counter] - tip_def_array[counter-1])/tip_def_array[counter]
    file_path_delete_neg = r"C:\Users\joepe\Documents\ProjectAero\OpenVSP-3.41.1-win64-Python3.11\OpenVSP-3.41.1-win64\python\update_geom_neg_DegenGeom.slc"
    if os.path.exists(file_path_delete_neg):
        os.remove(file_path_delete_neg)
        print("File deleted successfully.")
    else:
        print("The file does not exist.")
    # Call update geometry function
    update_geom_neg(wing_id)
    vsp.WriteVSPFile("update_geom_neg.vsp3")
    print("Twist and DeltaY updated successfully.")
    print("Negative load model saved successfully as update_geom_neg.vsp3")
    # Call VLM Function     # Creates pressure distribution file for the wing
    #vlm_fsi_neg()
    counter += 1
    print (counter)
    print (tip_deflection_error)

print ("counter",counter)
print ("tip def error",tip_deflection_error)
print ("tip def", tip_deflection)
#print(Ve_pos_init,alpha_e_pos_init,Ve_pos,alpha_e_pos)
"""
#LOAD DEFLECTED MODELS INTO THE CODE
"""

#update_geom_pos(wing_id)
#vsp.WriteVSPFile("update_geom_pos.vsp3")
#print("Twist and DeltaY updated successfully.")
#print("Positive load model saved successfully as update_geom_pos.vsp3")
#update_geom_neg(wing_id)
#vsp.WriteVSPFile("update_geom_neg.vsp3")
#print("Twist and DeltaY updated successfully.")
#print("Negative load model saved successfully as update_geom_neg.vsp3") 


#RUN VLM ON THE DEFLECTED MODELS

# Positive load vlm
#vlm_fsi_pos()

# Negative load vlm
#vlm_fsi_neg()

# Landing gear sizing – uses MTOW estimate for initial sizing 
#landing_gear_mass_array = landing_gear_sizing(ground_clearance, MTOW, nlg_distance_from_cg, mlg_distance_from_cg, sink_speed, z_cg)
#negative_fsi()
#positive_fsi()

# Write mass properties csv 
#write_aircraft_data_to_csv('aircraft_masses.csv',cabin_mass_array, wing_mass_array,vstab_mass_array, landing_gear_mass_array)

#vstab_mass_array,
