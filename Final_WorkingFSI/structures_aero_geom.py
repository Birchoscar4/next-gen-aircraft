import csv
import openvsp as vsp
from Initial_Geom import build_initial_geometry
from cabin_mass import compute_cabin_mass
from wing_model import compute_wing_mass
#from vstab_model import compute_vstab_mass
from landing_gear_sizing import landing_gear_sizing
from mass_properties import write_aircraft_data_to_csv
from structure_neg_fsi import update_geom_neg # Geometry Function from FSI
from structures_pos_fsi import update_geom_pos # Geometry Function from FSI
from vlm_fsi_neg import vlm_fsi_neg
from vlm_fsi_pos import vlm_fsi_pos

# Computes the generation of the initial geometry and the VSPAero analysis of the geometry
build_initial_geometry()

max_q_density = 0.35 #temp
max_q_velocity = 210 #temp
#Fuselage Inputs
frame_spacing = 0.2 #m
fuselage_sparcap_radius = 0.05 #m
#Wing Inputs
wing_sparcap_radius = 0.08
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

cabin_mass_array = compute_cabin_mass(frame_spacing, sparcap_radius=fuselage_sparcap_radius)

wing_mass_array = compute_wing_mass(spar_cap_radius=wing_sparcap_radius, num_stringers=wing_num_stringers, skin_thickness=wing_skin_thickness, spar_thickness=wing_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_wing, TE_spar=TE_spar_wing, file_path = "Prelim_2_DegenGeom.slc")

# Define Wing ID as the blended wing body aircraft wing
wing_id = vsp.FindGeomsWithName("WING")  # Replace "WING" with your actual component name
wing_id = wing_id[0]

update_geom_pos(wing_id)
vsp.WriteVSPFile("update_geom_pos.vsp3", vsp.SET_ALL)
print("Twist and DeltaY updated successfully.")
print("Positive load model saved successfully as update_geom_pos.vsp3")

# Save OpenVSP Just in case
vsp.Update()

# Run VLM code on the positive load case
vlm_fsi_pos()