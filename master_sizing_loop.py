#Aircraft Sizing Loop Master Code Draft 1
#Last update 02/04/2025 
#################LIST REQUIRED MODULES AND APIS HERE#######################
import pandas as pd
import math
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
pi = float(math.pi)
plt.ion()
import openvsp as vsp
from scipy.integrate import cumulative_trapezoid
import openmdao.api as om




#####################GEOMETRY AND CONFIGURATION CODE #####################




##########################AERODYNAMICS CODE #############################
#this is integrated into the fsi, it doesnt act as its own standalone chunk of code





##########################MISSION ANALYSIS CODE ############################# 
from openconcept.utilities import AddSubtractComp, DictIndepVarComp, plot_trajectory

# imports for the airplane model itself
from openconcept.aerodynamics import PolarDrag
from openconcept.examples.aircraft_data.BWB import data as acdata
from openconcept.mission import MissionWithReserve, IntegratorGroup
from openconcept.propulsion import RubberizedTurbofan
from openconcept.utilities import Integrator, AddSubtractComp, ElementMultiplyDivideComp, DictIndepVarComp
hydrogen = True

from BWB_Mission_With_Reserve.py import mission_plot
from BWB_Mission_With_Reserve.py import takeoff_and_landing

mission_plot = mission_plot()
takeoff_landing = takeoff_and_landing()






##########################STRUCTURES CODE #############################
import csv
from cabin_mass import compute_cabin_mass
from wing_model import compute_wing_mass
from vstab_model import compute_vstab_mass
from landing_gear_sizing import landing_gear_sizing
from mass_properties import write_aircraft_data_to_csv
from structure_neg_fsi import update_geom_neg # Geometry Function from FSI
from structures_pos_fsi import update_geom_pos # Geometry Function from FSI

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

# Fuselage curve fit equations, returns mass properties 
cabin_mass_array = compute_cabin_mass(frame_spacing, sparcap_radius=fuselage_sparcap_radius)
# aft_body_mass_array = compute_aft_body_mass()

# Vertical Stabiliser bending analysis
vstab_mass_array = compute_vstab_mass(spar_cap_radius=vstab_sparcap_radius, num_stringers=vstab_num_stringers, skin_thickness=vstab_skin_thickness, spar_thickness=vstab_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_vstab, TE_spar=TE_spar_vstab)

### FSI LOOPS HERE
tip_def_array = []
counter = 0

tip_deflection_error = 1
while tip_deflection_error > 0.05:
    # Call VLM Function     # Creates pressure distribution file for the wing
    
    # Call wing deflection function to produce deflection and twist distribution
    tip_deflection = compute_wing_mass(spar_cap_radius=wing_sparcap_radius, num_stringers=wing_num_stringers, skin_thickness=wing_skin_thickness, spar_thickness=wing_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_wing, TE_spar=TE_spar_wing)[9]
    # tip_def_array[counter] = tip_deflection
    tip_deflection_error = (tip_def_array[counter] - tip_def_array[counter-1])-tip_def_array[counter]
    # Call update geometry function
    counter += 1
    wing_mass_array = 
    update_geom_neg()


wing_mass_array = compute_wing_mass(spar_cap_radius=wing_sparcap_radius, num_stringers=wing_num_stringers, skin_thickness=wing_skin_thickness, spar_thickness=wing_spar_thickness,max_q_velocity=max_q_velocity, max_q_density=max_q_density,load_factors=load_factor, LE_spar=LE_spar_wing, TE_spar=TE_spar_wing)

update_geom_pos() 

# Landing gear sizing – uses MTOW estimate for initial sizing 
landing_gear_mass_array = landing_gear_sizing(ground_clearance, MTOW, nlg_distance_from_cg, mlg_distance_from_cg, sink_speed, z_cg)

# Write mass properties csv 
write_aircraft_data_to_csv('aircraft_masses.csv',cabin_mass_array, wing_mass_array, vstab_mass_array, landing_gear_mass_array)




##########################PROP & SYSTEMS CODE #############################
from openconcept.examples.aircraft_data.specific_BWB_sizing import data as acdata # this uses the specific updated values from the aero team.
from openconcept.aerodynamics import PolarDrag, ParasiteDragCoefficient_JetTransport, CleanCLmax, FlapCLmax
from openconcept.propulsion import RubberizedTurbofan
from openconcept.geometry import CylinderSurfaceArea, WingMACTrapezoidal
from openconcept.stability import HStabVolumeCoefficientSizing, VStabVolumeCoefficientSizing
from openconcept.weights import JetTransportEmptyWeight
from openconcept.mission import FullMissionWithReserve
from openconcept.utilities import Integrator, AddSubtractComp, ElementMultiplyDivideComp, DictIndepVarComp

class B738AircraftModel(om.Group):
    """
    A Boeing 737-800 aircraft model group. Instead of using known weight
    and drag estimates of the existing airplane, this group uses empirical
    weight and drag buildups to enable the design of clean sheet aircraft.
    """

    def initialize(self):
        self.options.declare("num_nodes", default=1, types=int, desc="Number of analysis points to run")
        self.options.declare("flight_phase", default=None, types=str, desc="Phase of mission this group lives in")

    def setup(self):
        nn = self.options["num_nodes"]
        phase = self.options["flight_phase"]
        in_takeoff = phase in ["v0v1", "v1v0", "v1vr", "rotate"]

        # ==============================================================================
        # Aerodynamics
        # ==============================================================================
        # -------------- Zero-lift drag coefficient buildup --------------
        drag_buildup_promotes = [
            "fltcond|Utrue",
            "fltcond|rho",
            "fltcond|T",
            "ac|geom|fuselage|length",
            "ac|geom|fuselage|height",
            "ac|geom|fuselage|S_wet",
            "ac|geom|hstab|S_ref",
            "ac|geom|hstab|AR",
            "ac|geom|hstab|taper",
            "ac|geom|hstab|toverc",
            "ac|geom|vstab|S_ref",
            "ac|geom|vstab|AR",
            "ac|geom|vstab|taper",
            "ac|geom|vstab|toverc",
            "ac|geom|wing|S_ref",
            "ac|geom|wing|AR",
            "ac|geom|wing|taper",
            "ac|geom|wing|toverc",
            "ac|geom|nacelle|length",
            "ac|geom|nacelle|S_wet",
            "ac|propulsion|num_engines",
        ]
        if in_takeoff:
            drag_buildup_promotes += ["ac|aero|takeoff_flap_deg", "ac|geom|wing|c4sweep"]
        self.add_subsystem(
            "zero_lift_drag",
            ParasiteDragCoefficient_JetTransport(num_nodes=nn, configuration="takeoff" if in_takeoff else "clean"),
            promotes_inputs=drag_buildup_promotes,
        )

        # -------------- Drag polar --------------
        self.add_subsystem(
            "drag_polar",
            PolarDrag(num_nodes=nn, vec_CD0=True),
            promotes_inputs=[
                "fltcond|CL",
                "fltcond|q",
                "ac|geom|wing|S_ref",
                "ac|geom|wing|AR",
                ("e", "ac|aero|polar|e"),
            ],
            promotes_outputs=["drag"],
        )
        self.connect("zero_lift_drag.CD0", "drag_polar.CD0")

        # ==============================================================================
        # Propulsion
        # ==============================================================================
        # -------------- N3 engine surrogate model --------------
        self.add_subsystem(
            "N3",
            RubberizedTurbofan(num_nodes=nn, engine="N3"),
            promotes_inputs=["throttle", "fltcond|h", "fltcond|M", "ac|propulsion|engine|rating"],
        )

        # -------------- Multiply fuel flow and thrust by the number of active engines --------------
        # propulsor_active is 0 if failed engine and 1 otherwise, so
        # num active engines = num engines - 1 + propulsor_active
        self.add_subsystem(
            "num_engine_calc",
            AddSubtractComp(
                output_name="num_active_engines",
                input_names=["num_engines", "propulsor_active", "one"],
                vec_size=[1, nn, 1],
                scaling_factors=[1, 1, -1],
            ),
            promotes_inputs=[("num_engines", "ac|propulsion|num_engines"), "propulsor_active"],
        )
        self.set_input_defaults("num_engine_calc.one", 1.0)

        prop_mult = self.add_subsystem(
            "propulsion_multiplier", ElementMultiplyDivideComp(), promotes_outputs=["thrust"]
        )
        prop_mult.add_equation(
            output_name="thrust",
            input_names=["thrust_per_engine", "num_active_engines_1"],
            vec_size=nn,
            input_units=["lbf", None],
        )
        prop_mult.add_equation(
            output_name="fuel_flow",
            input_names=["fuel_flow_per_engine", "num_active_engines_2"],
            vec_size=nn,
            input_units=["kg/s", None],
        )
        self.connect("N3.fuel_flow", "propulsion_multiplier.fuel_flow_per_engine")
        self.connect("N3.thrust", "propulsion_multiplier.thrust_per_engine")

        # This hacky thing is necessary to enable two equations to pull from the same input
        self.connect(
            "num_engine_calc.num_active_engines",
            ["propulsion_multiplier.num_active_engines_1", "propulsion_multiplier.num_active_engines_2"],
        )

        # ==============================================================================
        # Weight
        # ==============================================================================
        # -------------- Integrate fuel burn --------------
        integ = self.add_subsystem(
            "fuel_burn_integ", Integrator(num_nodes=nn, diff_units="s", method="simpson", time_setup="duration")
        )
        integ.add_integrand(
            "fuel_burn",
            rate_name="fuel_flow",
            rate_units="kg/s",
            lower=0.0,
            upper=1e6,
        )
        self.connect("propulsion_multiplier.fuel_flow", "fuel_burn_integ.fuel_flow")

        # -------------- Subtract fuel burn from takeoff weight --------------
        self.add_subsystem(
            "weight_calc",
            AddSubtractComp(
                output_name="weight",
                input_names=["ac|weights|MTOW", "fuel_burn"],
                scaling_factors=[1, -1],
                vec_size=[1, nn],
                units="kg",
            ),
            promotes_inputs=["ac|weights|MTOW"],
            promotes_outputs=["weight"],
        )
        self.connect("fuel_burn_integ.fuel_burn", "weight_calc.fuel_burn")

class B738SizingMissionAnalysis(om.Group):
    """
    Group that performs mission analysis, geometry calculations, and operating empty weight estimate.
    """

    def initialize(self):
        self.options.declare("num_nodes", default=11, types=int, desc="Analysis points per mission phase")

    def setup(self):
        nn = self.options["num_nodes"]

        # ==============================================================================
        # Variables from B738_sizing data file
        # ==============================================================================
        dv = self.add_subsystem("ac_vars", DictIndepVarComp(acdata), promotes_outputs=["*"])
        dv_outputs = [
            # -------------- Aero --------------
            "ac|aero|polar|e",
            "ac|aero|Mach_max",
            "ac|aero|Vstall_land",
            "ac|aero|airfoil_Cl_max",
            "ac|aero|takeoff_flap_deg",
            # -------------- Propulsion --------------
            "ac|propulsion|engine|rating",
            "ac|propulsion|num_engines",
            # -------------- Geometry --------------
            # Wing
            "ac|geom|wing|S_ref",
            "ac|geom|wing|AR",
            "ac|geom|wing|c4sweep",
            "ac|geom|wing|taper",
            "ac|geom|wing|toverc",
            # Horizontal stabilizer
            "ac|geom|hstab|AR",
            "ac|geom|hstab|c4sweep",
            "ac|geom|hstab|taper",
            "ac|geom|hstab|toverc",
            # Vertical stabilizer
            "ac|geom|vstab|AR",
            "ac|geom|vstab|c4sweep",
            "ac|geom|vstab|taper",
            "ac|geom|vstab|toverc",
            # Fuselage
            "ac|geom|fuselage|length",
            "ac|geom|fuselage|height",
            # Nacelle
            "ac|geom|nacelle|length",
            "ac|geom|nacelle|diameter",
            # Main gear
            "ac|geom|maingear|length",
            "ac|geom|maingear|num_wheels",
            "ac|geom|maingear|num_shock_struts",
            # Nose gear
            "ac|geom|nosegear|length",
            "ac|geom|nosegear|num_wheels",
            # -------------- Weights --------------
            "ac|weights|W_payload",
            # -------------- Miscellaneous --------------
            "ac|num_passengers_max",
            "ac|num_flight_deck_crew",
            "ac|num_cabin_crew",
            "ac|cabin_pressure",
        ]
        for output_name in dv_outputs:
            dv.add_output_from_dict(output_name)

        # ==============================================================================
        # Geometry
        # ==============================================================================
        # -------------- Estimate wing to tail quarter chord as half fuselage length --------------
        self.add_subsystem(
            "tail_lever_arm_estimate",
            AddSubtractComp(
                output_name="c4_to_wing_c4",
                input_names=["fuselage_length"],
                units="m",
                scaling_factors=[0.5],
            ),
            promotes_inputs=[("fuselage_length", "ac|geom|fuselage|length")],
        )
        self.connect(
            "tail_lever_arm_estimate.c4_to_wing_c4", ["ac|geom|hstab|c4_to_wing_c4", "ac|geom|vstab|c4_to_wing_c4"]
        )

        # -------------- Compute mean aerodynamic chord assuming trapezoidal wing --------------
        self.add_subsystem(
            "wing_MAC",
            WingMACTrapezoidal(),
            promotes_inputs=[
                ("S_ref", "ac|geom|wing|S_ref"),
                ("AR", "ac|geom|wing|AR"),
                ("taper", "ac|geom|wing|taper"),
            ],
            promotes_outputs=[("MAC", "ac|geom|wing|MAC")],
        )

        # -------------- Vertical and horizontal tail area --------------
        self.add_subsystem(
            "vstab_area",
            VStabVolumeCoefficientSizing(),
            promotes_inputs=["ac|geom|wing|S_ref", "ac|geom|wing|AR", "ac|geom|vstab|c4_to_wing_c4"],
            promotes_outputs=["ac|geom|vstab|S_ref"],
        )
        self.add_subsystem(
            "hstab_area",
            HStabVolumeCoefficientSizing(),
            promotes_inputs=["ac|geom|wing|S_ref", "ac|geom|wing|MAC", "ac|geom|hstab|c4_to_wing_c4"],
            promotes_outputs=["ac|geom|hstab|S_ref"],
        )

        # -------------- Compute the fuselage and nacelle wetted areas assuming a cylinder --------------
        self.add_subsystem(
            "nacelle_wetted_area",
            CylinderSurfaceArea(),
            promotes_inputs=[("L", "ac|geom|nacelle|length"), ("D", "ac|geom|nacelle|diameter")],
            promotes_outputs=[("A", "ac|geom|nacelle|S_wet")],
        )
        self.add_subsystem(
            "fuselage_wetted_area",
            CylinderSurfaceArea(),
            promotes_inputs=[("L", "ac|geom|fuselage|length"), ("D", "ac|geom|fuselage|height")],
            promotes_outputs=[("A", "ac|geom|fuselage|S_wet")],
        )

        # ==============================================================================
        # Operating empty weight
        # ==============================================================================
        # Estimate MLW as 80% of MTOW so it's not necessary to know apriori.
        # MLW is used only for the landing gear weight estimate, so it's accuracy isn't hugely important.
        self.add_subsystem(
            "MLW_calc",
            AddSubtractComp(
                output_name="ac|weights|MLW",
                input_names=["ac|weights|MTOW"],
                units="kg",
                scaling_factors=[0.8],
            ),
            promotes_inputs=["ac|weights|MTOW"],
            promotes_outputs=["ac|weights|MLW"],
        )

        self.add_subsystem(
            "empty_weight",
            JetTransportEmptyWeight(),
            promotes_inputs=[
                "ac|num_passengers_max",
                "ac|num_flight_deck_crew",
                "ac|num_cabin_crew",
                "ac|cabin_pressure",
                "ac|aero|Mach_max",
                "ac|aero|Vstall_land",
                "ac|geom|wing|S_ref",
                "ac|geom|wing|AR",
                "ac|geom|wing|c4sweep",
                "ac|geom|wing|taper",
                "ac|geom|wing|toverc",
                "ac|geom|hstab|S_ref",
                "ac|geom|hstab|AR",
                "ac|geom|hstab|c4sweep",
                "ac|geom|hstab|c4_to_wing_c4",
                "ac|geom|vstab|S_ref",
                "ac|geom|vstab|AR",
                "ac|geom|vstab|c4sweep",
                "ac|geom|vstab|toverc",
                "ac|geom|vstab|c4_to_wing_c4",
                "ac|geom|fuselage|height",
                "ac|geom|fuselage|length",
                "ac|geom|fuselage|S_wet",
                "ac|geom|maingear|length",
                "ac|geom|maingear|num_wheels",
                "ac|geom|maingear|num_shock_struts",
                "ac|geom|nosegear|length",
                "ac|geom|nosegear|num_wheels",
                "ac|propulsion|engine|rating",
                "ac|propulsion|num_engines",
                "ac|weights|MTOW",
                "ac|weights|MLW",
            ],
            promotes_outputs=[("OEW", "ac|weights|OEW")],
        )

        # ==============================================================================
        # CL max in cruise and takeoff
        # ==============================================================================
        self.add_subsystem(
            "CL_max_cruise",
            CleanCLmax(),
            promotes_inputs=["ac|aero|airfoil_Cl_max", "ac|geom|wing|c4sweep"],
            promotes_outputs=[("CL_max_clean", "ac|aero|CLmax_cruise")],
        )
        self.add_subsystem(
            "CL_max_takeoff",
            FlapCLmax(),
            promotes_inputs=[
                ("flap_extension", "ac|aero|takeoff_flap_deg"),
                "ac|geom|wing|c4sweep",
                "ac|geom|wing|toverc",
                ("CL_max_clean", "ac|aero|CLmax_cruise"),
            ],
            promotes_outputs=[("CL_max_flap", "ac|aero|CLmax_TO")],
        )

        # ==============================================================================
        # Remaining misc stuff and input settings
        # ==============================================================================
        # -------------- Set MTOW to OEW + payload + fuel burn --------------
        self.add_subsystem(
            "MTOW_calc",
            AddSubtractComp(
                output_name="MTOW",
                input_names=["OEW", "W_payload", "W_fuel"],
                units="kg",
                lower=1e-6,
            ),
            promotes_inputs=[("OEW", "ac|weights|OEW"), ("W_payload", "ac|weights|W_payload")],
            promotes_outputs=[("MTOW", "ac|weights|MTOW")],
        )
        self.connect(
            "mission.loiter.fuel_burn_integ.fuel_burn_final",
            ["MTOW_calc.W_fuel", "empty_weight.ac|weights|W_fuel_max"],
        )

        # -------------- Initial guesses for important solver states for better performance --------------
        self.set_input_defaults("ac|weights|MTOW", 50e3, units="kg")
        for fuel_burn_var in ["MTOW_calc.W_fuel", "empty_weight.ac|weights|W_fuel_max"]:
            self.set_input_defaults(fuel_burn_var, 30e3, units="kg")

        # ==============================================================================
        # Mission analysis
        # ==============================================================================
        self.add_subsystem(
            "mission",
            FullMissionWithReserve(num_nodes=nn, aircraft_model=B738AircraftModel),
            promotes_inputs=["ac|*"],
        )

def set_mission_profile(prob):
    """
    Set the parameters in the OpenMDAO problem that define the mission profile.

    Parameters
    ----------
    prob : OpenMDAO Problem
        Problem with B378MissionAnalysis as model in which to set values
    """
    # Get the number of nodes in the mission
    nn = prob.model.options["num_nodes"]

    # ==============================================================================
    # Basic mission phases
    # ==============================================================================
    # -------------- Climb --------------
    prob.set_val("mission.climb.fltcond|vs", np.linspace(2300.0, 400.0, nn), units="ft/min")
    prob.set_val("mission.climb.fltcond|Ueas", np.linspace(230, 252, nn), units="kn")

    # -------------- Cruise --------------
    prob.set_val("mission.cruise.fltcond|vs", np.full((nn,), 0.0), units="ft/min")
    prob.set_val("mission.cruise.fltcond|Ueas", np.linspace(252, 252, nn), units="kn")

    # -------------- Descent --------------
    prob.set_val("mission.descent.fltcond|vs", np.linspace(-1300, -800, nn), units="ft/min")
    prob.set_val("mission.descent.fltcond|Ueas", np.linspace(252, 250, nn), units="kn")

    # ==============================================================================
    # Reserve mission phases
    # ==============================================================================
    # -------------- Reserve climb --------------
    prob.set_val("mission.reserve_climb.fltcond|vs", np.linspace(3000.0, 2300.0, nn), units="ft/min")
    prob.set_val("mission.reserve_climb.fltcond|Ueas", np.linspace(230, 230, nn), units="kn")

    # -------------- Reserve cruise --------------
    prob.set_val("mission.reserve_cruise.fltcond|vs", np.full((nn,), 0.0), units="ft/min")
    prob.set_val("mission.reserve_cruise.fltcond|Ueas", np.linspace(250, 250, nn), units="kn")

    # -------------- Reserve descent --------------
    prob.set_val("mission.reserve_descent.fltcond|vs", np.linspace(-800, -800, nn), units="ft/min")
    prob.set_val("mission.reserve_descent.fltcond|Ueas", np.full((nn,), 250.0), units="kn")

    # -------------- Loiter --------------
    prob.set_val("mission.loiter.fltcond|vs", np.linspace(0.0, 0.0, nn), units="ft/min")
    prob.set_val("mission.loiter.fltcond|Ueas", np.full((nn,), 250.0), units="kn")

    # ==============================================================================
    # Other parameters
    # ==============================================================================
    prob.set_val("mission.cruise|h0", 35000.0, units="ft")
    prob.set_val("mission.reserve|h0", 15000.0, units="ft")
    prob.set_val("mission.mission_range", 1600, units="nmi") #edited to consider the 1600nmi mission specification range

    # -------------- Set takeoff speed guesses to improve solver performance --------------
    prob.set_val("mission.v0v1.fltcond|Utrue", np.full((nn,), 100.0), units="kn")
    prob.set_val("mission.v1vr.fltcond|Utrue", np.full((nn,), 100.0), units="kn")
    prob.set_val("mission.v1v0.fltcond|Utrue", np.full((nn,), 100.0), units="kn")

    # Converge the model first with an easier mission profile and work up to the intended
    # mission profile. This is needed to help the Newton solver converge the actual mission.
    prob.set_val("mission.descent.fltcond|vs", np.linspace(-800, -800, nn), units="ft/min")
    prob.set_val("mission.cruise|h0", 5000.0, units="ft")
    prob.set_val("mission.reserve|h0", 1000.0, units="ft")
    prob.set_val("mission.mission_range", 500, units="nmi")
    prob.set_val("mission.reserve_range", 100, units="nmi")
    prob.run_model()

    # Almost there, just not quite with descent rate
    prob.set_val("mission.cruise|h0", 35000.0, units="ft")
    prob.set_val("mission.reserve|h0", 15000.0, units="ft")
    prob.set_val("mission.mission_range", 1600, units="nmi") #edited to consider the 1600nmi mission specification range
    prob.set_val("mission.reserve_range", 200, units="nmi") #edited from 200 to comply with the ICAO minimum fuel requirement.
    prob.run_model()

    # Finally, set the descent rate we want
    prob.set_val("mission.descent.fltcond|vs", np.linspace(-1300, -800, nn), units="ft/min")

def plot_results(prob, filename=None):
    """
    Make a plot with the results of the mission analysis.

    Parameters
    ----------
    prob : OpenMDAO Problem
        Problem with B738SizingMissionAnalysis model that has been run
    filename : str (optional)
        Filename to save to, by default will show plot
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    axs = axs.flatten()

    for phase in ["climb", "cruise", "descent", "reserve_climb", "reserve_cruise", "reserve_descent", "loiter"]:
        dist = prob.get_val(f"mission.{phase}.range", units="nmi")

        axs[0].plot(dist, prob.get_val(f"mission.{phase}.fltcond|h", units="ft"), color="tab:red")
        axs[1].plot(dist, prob.get_val(f"mission.{phase}.fltcond|M"), color="tab:blue")
        axs[2].plot(dist, prob.get_val(f"mission.{phase}.fltcond|vs", units="ft/min"), color="tab:blue")
        axs[3].plot(dist, prob.get_val(f"mission.{phase}.weight", units="kg"), color="tab:blue")
        axs[4].plot(dist, prob.get_val(f"mission.{phase}.drag", units="lbf"), color="tab:blue")
        axs[4].plot(dist, prob.get_val(f"mission.{phase}.thrust", units="lbf"), color="tab:orange")
        axs[5].plot(dist, prob.get_val(f"mission.{phase}.throttle") * 100, color="tab:blue")

    axs[0].set_ylabel("Altitude (ft)")
    axs[1].set_ylabel("Mach number")
    axs[2].set_ylabel("Vertical speed (ft/min)")
    axs[3].set_ylabel("Mass (kg)")
    axs[4].set_ylabel("Longitudinal force (lb)")
    axs[5].set_ylabel("Throttle (%)")
    axs[4].legend(["Drag", "Thrust"])

    #Extract the thrust, Mach number and range2 values and output them as an excel file for all mission phases.
    thrust_climb = prob.get_val(f"mission.{"climb"}.thrust")
    thrust_cruise = prob.get_val(f"mission.{"cruise"}.thrust")
    thrust_descent = prob.get_val(f"mission.{"descent"}.thrust")    
    thrust_reserve_climb = prob.get_val(f"mission.{"reserve_climb"}.thrust")
    thrust_reserve_cruise = prob.get_val(f"mission.{"reserve_cruise"}.thrust")
    thrust_reserve_descent = prob.get_val(f"mission.{"reserve_descent"}.thrust")
    thrust_loiter = prob.get_val(f"mission.{"loiter"}.thrust")
    df = pd.DataFrame(thrust_climb)
    df.to_excel('thrust_climb.xlsx', index=False)
    df = pd.DataFrame(thrust_cruise)
    df.to_excel('thrust_cruise.xlsx', index=False)
    df = pd.DataFrame(thrust_descent)
    df.to_excel('thrust_descent.xlsx', index=False)
    df = pd.DataFrame(thrust_reserve_climb)
    df.to_excel('thrust_reserve_climb.xlsx', index=False)
    df = pd.DataFrame(thrust_reserve_cruise)
    df.to_excel('thrust_reserve_cruise.xlsx', index=False)
    df = pd.DataFrame(thrust_reserve_descent)
    df.to_excel('thrust_reserve_descent.xlsx', index=False)
    df = pd.DataFrame(thrust_loiter)
    df.to_excel('thrust_loiter.xlsx', index=False)

    mach_no_climb = prob.get_val(f"mission.{"climb"}.fltcond|M")
    mach_no_cruise = prob.get_val(f"mission.{"cruise"}.fltcond|M")
    mach_no_descent = prob.get_val(f"mission.{"descent"}.fltcond|M")    
    mach_no_reserve_climb = prob.get_val(f"mission.{"reserve_climb"}.fltcond|M")
    mach_no_reserve_cruise = prob.get_val(f"mission.{"reserve_cruise"}.fltcond|M")
    mach_no_reserve_descent = prob.get_val(f"mission.{"reserve_descent"}.fltcond|M")
    mach_no_loiter = prob.get_val(f"mission.{"loiter"}.fltcond|M")
    df = pd.DataFrame(mach_no_climb)
    df.to_excel('mach_no_climb.xlsx', index=False)
    df = pd.DataFrame(mach_no_cruise)
    df.to_excel('mach_no_cruise.xlsx', index=False)
    df = pd.DataFrame(mach_no_descent)
    df.to_excel('mach_no_descent.xlsx', index=False)
    df = pd.DataFrame(mach_no_reserve_climb)
    df.to_excel('mach_no_reserve_climb.xlsx', index=False)
    df = pd.DataFrame(mach_no_reserve_cruise)
    df.to_excel('mach_no_reserve_cruise.xlsx', index=False)
    df = pd.DataFrame(mach_no_reserve_descent)
    df.to_excel('mach_no_reserve_descent.xlsx', index=False)
    df = pd.DataFrame(mach_no_loiter)
    df.to_excel('mach_no_loiter.xlsx', index=False)

    range2_climb = prob.get_val(f"mission.{"climb"}.range")
    range2_cruise = prob.get_val(f"mission.{"cruise"}.range")
    range2_descent = prob.get_val(f"mission.{"descent"}.range")    
    range2_reserve_climb = prob.get_val(f"mission.{"reserve_climb"}.range")
    range2_reserve_cruise = prob.get_val(f"mission.{"reserve_cruise"}.range")
    range2_reserve_descent = prob.get_val(f"mission.{"reserve_descent"}.range")
    range2_loiter = prob.get_val(f"mission.{"loiter"}.range")
    df = pd.DataFrame(range2_climb)
    df.to_excel('range2_climb.xlsx', index=False)
    df = pd.DataFrame(range2_cruise)
    df.to_excel('range2_cruise.xlsx', index=False)
    df = pd.DataFrame(range2_descent)
    df.to_excel('range2_descent.xlsx', index=False)
    df = pd.DataFrame(range2_reserve_climb)
    df.to_excel('range2_reserve_climb.xlsx', index=False)
    df = pd.DataFrame(range2_reserve_cruise)
    df.to_excel('range2_reserve_cruise.xlsx', index=False)
    df = pd.DataFrame(range2_reserve_descent)
    df.to_excel('range2_reserve_descent.xlsx', index=False)
    df = pd.DataFrame(range2_loiter)
    df.to_excel('range2_loiter.xlsx', index=False)

    for i in range(6):
        axs[i].set_xlabel("Distance flown (nmi)")
        axs[i].spines[["right", "top"]].set_visible(False)
        if i != 1:
            axs[i].get_yaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ",")))
        axs[i].get_xaxis().set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ",")))

    plt.tight_layout()

    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)

def run_738_sizing_analysis(num_nodes=21):
    p = om.Problem()
    p.model = B738SizingMissionAnalysis(num_nodes=num_nodes)

    # -------------- Add solvers --------------
    p.model.nonlinear_solver = om.NewtonSolver()
    p.model.nonlinear_solver.options["iprint"] = 2
    p.model.nonlinear_solver.options["solve_subsystems"] = True
    p.model.nonlinear_solver.options["maxiter"] = 20
    p.model.nonlinear_solver.options["atol"] = 1e-9
    p.model.nonlinear_solver.options["rtol"] = 1e-9

    p.model.linear_solver = om.DirectSolver()

    p.model.nonlinear_solver.linesearch = om.BoundsEnforceLS()
    p.model.nonlinear_solver.linesearch.options["print_bound_enforce"] = False

    p.setup()

    set_mission_profile(p)

    p.run_model()
    return p

if __name__ == "__main__":
    p = run_738_sizing_analysis()
    om.n2(p, show_browser=False)

    # Print some useful numbers
    print("\n\n================= Computed values =================")
    print(f"MTOW: {p.get_val('ac|weights|MTOW', units='lb').item():.1f} lb")
    print(f"Payload weight: {p.get_val('ac|weights|W_payload', units='lb').item()} lb")
    print(f"OEW: {p.get_val('ac|weights|OEW', units='lb').item():.1f} lb")
    print(f"Fuel burned: {p.get_val('mission.descent.fuel_burn_integ.fuel_burn_final', units='lb').item():.1f} lb")
    print(f"CL max cruise: {p.get_val('ac|aero|CLmax_cruise').item():.3f}")
    print(f"CL max takeoff: {p.get_val('ac|aero|CLmax_TO').item():.3f}")
    print(f"Balanced field length (continue): {p.get_val('mission.bfl.distance_continue', units='ft').item():.1f} ft")
    print(
        f"Balanced field length (abort): {p.get_val('mission.bfl.distance_abort', units='ft').item():.1f} ft (this should be the same as continue)"
    )

print("Thrust has been written to 'thrust_BWB.xlsx' and Mach number has been written to 'MachNo_BWB.xlsx.")

plot_results(p, filename="BWB_mission_profile.pdf")

## LH2 mass to volume conversion code
# User inputs
radius_constraint = 1.5 #This is the maximum radius of the tank due to fuselage available space in m
tank_quantity = 2 #This determines what the overall tank volume is divided by
MLI_thickness = 0.1016 #This is the insulation thickness in m, it is inputted in to the code after a conversion to inches

LH2_mass_kg = p.get_val('mission.descent.fuel_burn_integ.fuel_burn_final', units='lb') * 0.453592 #Convert LH2 mass in lb to kg
print("Mass of LH2 required:", LH2_mass_kg, "kg")
LH2_density = 70.9 #In kg/m^3
LH2_volume = LH2_mass_kg/LH2_density #Finds volume of LH2 required

LH2_tank_volume = LH2_volume/0.93 #This accounts for the 93% maximum fill rate of LH2 in a tank in m^3
tank_length = LH2_tank_volume/(radius_constraint**2*np.pi)/tank_quantity #This calculates the corresponding length of the tank subject to the required volume, radius constraint and tank quantity

# Outputs for the tank sizing code to use, conversions to the correct units
radius_constraint_ft = radius_constraint*3.28084
tank_length_ft = tank_length*3.28084
MLI_thickness_inch = MLI_thickness*39.3701

print("Tank radius:", radius_constraint_ft, "ft")
print("Tank length:", tank_length_ft, "ft")
print("MLI thickness:", MLI_thickness_inch, "in")

## Tank weight code
class VacuumTankWeight(om.Group):
    """
    Sizes the structure and computes the weight of the tank's vacuum walls.
    This includes the weight of MLI.

    .. code-block:: text

              |--- length ---|
             . -------------- .         ---
          ,'                    `.       | radius
         /                        \      |
        |                          |    ---
         \                        /
          `.                    ,'
             ` -------------- '

    Inputs
    ------
    environment_design_pressure : float
        Maximum environment exterior pressure expected, probably ~1 atmosphere (scalar, Pa)
    max_expected_operating_pressure : float
        Maximum expected operating pressure of tank (scalar, Pa)
    vacuum_gap : float
        Thickness of vacuum gap, used to compute radius of outer vacuum wall (scalar, m)
    radius : float
        Tank inner radius of the cylinder and hemispherical end caps (scalar, m)
    length : float
        Length of JUST THE CYLIDRICAL part of the tank (scalar, m)
    N_layers : float
        Number of reflective sheild layers in the MLI, should be at least ~10 for model
        to retain reasonable accuracy (scalar, dimensionless)

    Outputs
    -------
    weight : float
        Weight of the tank walls (scalar, kg)

    Options
    -------
    weight_fudge_factor : float
        Multiplier on tank weight to account for supports, valves, etc., by default 1.1
    stiffening_multiplier : float
        Machining stiffeners into the inner side of the vacuum shell enhances its buckling
        performance, enabling weight reductions. The value provided in this option is a
        multiplier on the outer wall thickness. The default value of 0.8 is higher than it
        would be if it were purely empirically determined from Sullivan et al. 2006
        (https://ntrs.nasa.gov/citations/20060021606), but has been made much more
        conservative to fall more in line with ~60% gravimetric efficiency tanks
    inner_safety_factor : float
        Safety factor for sizing inner wall, by default 1.5
    inner_yield_stress : float
        Yield stress of inner wall material (Pa), by default Al 2014-T6 taken from Table IV of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    inner_density : float
        Density of inner wall material (kg/m^3), by default Al 2014-T6 taken from Table IV of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    outer_safety_factor : float
        Safety factor for sizing outer wall, by default 2
    outer_youngs_modulus : float
        Young's modulus of outer wall material (Pa), by default LiAl 2090 taken from Table XIII of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    outer_density : float
        Density of outer wall material (kg/m^3), by default LiAl 2090 taken from Table XIII of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    """

    def initialize(self):
        self.options.declare("weight_fudge_factor", default=1.1, desc="Weight multiplier to account for other stuff")
        self.options.declare("stiffening_multiplier", default=0.8, desc="Multiplier on wall thickness")
        self.options.declare("inner_safety_factor", default=1.5, desc="Safety factor on inner wall thickness")
        self.options.declare("inner_yield_stress", default=413.7e6, desc="Yield stress of inner wall material in Pa")
        self.options.declare("inner_density", default=2796.0, desc="Density of inner wall material in kg/m^3")
        self.options.declare("outer_safety_factor", default=2.0, desc="Safety factor on outer wall thickness")
        self.options.declare("outer_youngs_modulus", default=8.0e10, desc="Young's modulus of outer wall material, Pa")
        self.options.declare("outer_density", default=2699.0, desc="Density of outer wall material in kg/m^3")

    def setup(self):
        # Inner tank wall thickness and weight computation
        self.add_subsystem(
            "inner_wall",
            PressureVesselWallThickness(
                safety_factor=self.options["inner_safety_factor"],
                yield_stress=self.options["inner_yield_stress"],
                density=self.options["inner_density"],
            ),
            promotes_inputs=[("design_pressure_differential", "max_expected_operating_pressure"), "radius", "length"],
        )

        # Compute radius of outer tank wall
        self.add_subsystem(
            "outer_radius",
            AddSubtractComp(
                output_name="outer_radius",
                input_names=["radius", "vacuum_gap"],
                scaling_factors=[1, 1],
                lower=0.0,
                units="m",
            ),
            promotes_inputs=["radius", "vacuum_gap"],
        )

        # Outer tank wall thickness and weight computation
        self.add_subsystem(
            "outer_wall",
            VacuumWallThickness(
                safety_factor=self.options["outer_safety_factor"],
                stiffening_multiplier=self.options["stiffening_multiplier"],
                youngs_modulus=self.options["outer_youngs_modulus"],
                density=self.options["outer_density"],
            ),
            promotes_inputs=[("design_pressure_differential", "environment_design_pressure"), "length"],
        )
        self.connect("outer_radius.outer_radius", "outer_wall.radius")

        # Compute the weight of the MLI
        self.add_subsystem("MLI", MLIWeight(), promotes_inputs=["radius", "length", "N_layers"])

        # Compute total weight multiplied by fudge factor and tank quantity chosen
        W_mult = self.options["weight_fudge_factor"]
        self.add_subsystem(
            "total_weight",
            AddSubtractComp(
                output_name="weight",
                input_names=["W_outer", "W_inner", "W_MLI"],
                scaling_factors=[W_mult, W_mult, W_mult],
                lower=0.0,
                units="kg",
            ),
            promotes_outputs=["weight"]*tank_quantity, # Single tank mass multiplied by the tank quantity chosen
        )
        self.connect("inner_wall.weight", "total_weight.W_inner")
        self.connect("outer_wall.weight", "total_weight.W_outer")
        self.connect("MLI.weight", "total_weight.W_MLI")

        # Set defaults for inputs promoted from multiple sources
        self.set_input_defaults("radius", radius_constraint, units="m") #by default 1m, 1m also used for the BWB project analysis for cylindrical tank analysis scenario 1.
        self.set_input_defaults("length", tank_length, units="m") # by default 0.5m, this is the length of purely the cylindrical part of the tank. 

class PressureVesselWallThickness(om.ExplicitComponent):
    """
    Compute the wall thickness of a metallic pressure vessel to support a specified
    pressure load. The model assumes an isotropic wall material, hence the metallic
    constraint. This uses a simple equation to compute the hoop stress (also referred
    to as Barlow's formula) to size the wall thickness.

    This component assumes that the wall is thin enough relative to the radius such that
    it is valid to compute the weight as the product of the surface area, wall thickness,
    and material density.

    .. code-block:: text

              |--- length ---|
             . -------------- .         ---
          ,'                    `.       | radius
         /                        \      |
        |                          |    ---
         \                        /
          `.                    ,'
             ` -------------- '

    Inputs
    ------
    design_pressure_differential : float
        The maximum pressure differential between the interior and exterior of the
        pressure vessel that is used to size the wall thickness; should ALWAYS
        be positive, otherwise wall thickness and weight will be negative (scalar, Pa)
    radius : float
        Inner radius of the cylinder and hemispherical end caps (scalar, m)
    length : float
        Length of JUST THE CYLIDRICAL part of the tank (scalar, m)

    Outputs
    -------
    thickness : float
        Pressure vessel wall thickness (scalar, m)
    weight : float
        Weight of the wall (scalar, kg)

    Options
    -------
    safety_factor : float
        Safety factor for sizing wall, by default 2
    yield_stress : float
        Yield stress of wall material (Pa), by default LiAl 2090 taken from Table XIII of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    density : float
        Density of wall material (kg/m^3), by default LiAl 2090 taken from Table XIII of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    """

    def initialize(self):
        self.options.declare("safety_factor", default=2.0, desc="Safety factor on wall thickness") # conservative safety factor frequently used in the Aerospace industry.
        self.options.declare("yield_stress", default=470.2e6, desc="Yield stress of wall material in Pa") # yield stress of an aviation grade aluminium alloy.
        self.options.declare("density", default=2699.0, desc="Density of wall material in kg/m^3") # density of an aviation grade aluminium alloy.

    def setup(self):
        self.add_input("design_pressure_differential", val=3e5, units="Pa") # difference in pressure between the LH2 in the tank and the external pressure. By default 30,000Pa.
        self.add_input("radius", val=radius_constraint, units="m") # by default 0.5m, 1m used for this initial cylindrical tank analysis. 
        self.add_input("length", val=tank_length, units="m") #by default 2m, 10m used for this initial cylindrical tank analysis.

        self.add_output("thickness", lower=0.0, units="m")
        self.add_output("weight", lower=0.0, units="kg")

        self.declare_partials("thickness", ["design_pressure_differential", "radius"])
        self.declare_partials("weight", ["design_pressure_differential", "radius", "length"])

    def compute(self, inputs, outputs):
        p = inputs["design_pressure_differential"]
        r = inputs["radius"]
        L = inputs["length"]
        SF = self.options["safety_factor"]
        yield_stress = self.options["yield_stress"]
        density = self.options["density"]

        outputs["thickness"] = p * r * SF / yield_stress

        surface_area = 4 * np.pi * r**2 + 2 * np.pi * r * L
        outputs["weight"] = surface_area * outputs["thickness"] * density

    def compute_partials(self, inputs, J):
        p = inputs["design_pressure_differential"]
        r = inputs["radius"]
        L = inputs["length"]
        SF = self.options["safety_factor"]
        yield_stress = self.options["yield_stress"]
        density = self.options["density"]

        t = p * r * SF / yield_stress

        J["thickness", "design_pressure_differential"] = r * SF / yield_stress
        J["thickness", "radius"] = p * SF / yield_stress

        A = 4 * np.pi * r**2 + 2 * np.pi * r * L
        dAdr = 8 * np.pi * r + 2 * np.pi * L
        dAdL = 2 * np.pi * r
        J["weight", "design_pressure_differential"] = A * J["thickness", "design_pressure_differential"] * density
        J["weight", "radius"] = (dAdr * t + A * J["thickness", "radius"]) * density
        J["weight", "length"] = dAdL * t * density

class VacuumWallThickness(om.ExplicitComponent):
    """
    Compute the wall thickness when the exterior pressure is greater than the interior
    one. This applies to the outer wall of a vacuum-insulated tank. It does this by
    computing the necessary wall thickness for a cylindrical shell under uniform compression
    and sphere under uniform compression and taking the maximum thickness of the two.

    The equations are from Table 15.2 of Roark's Formulas for Stress and Strain, 9th
    Edition by Budynas and Sadegh.

    This component assumes that the wall is thin relative to the radius.

    .. code-block:: text

              |--- length ---|
             . -------------- .         ---
          ,'                    `.       | radius
         /                        \      |
        |                          |    ---
         \                        /
          `.                    ,'
             ` -------------- '

    Inputs
    ------
    design_pressure_differential : float
        The maximum pressure differential between the interior and exterior of the
        pressure vessel that is used to size the wall thickness; should ALWAYS
        be positive (scalar, Pa)
    radius : float
        Inner radius of the cylinder and hemispherical end caps (scalar, m)
    length : float
        Length of JUST THE CYLIDRICAL part of the tank (scalar, m)

    Outputs
    -------
    thickness : float
        Pressure vessel wall thickness (scalar, m)
    weight : float
        Weight of the wall (scalar, kg)

    Options
    -------
    safety_factor : float
        Safety factor for sizing wall applied to design pressure, by default 2
    stiffening_multiplier : float
        Machining stiffeners into the inner side of the vacuum shell enhances its buckling
        performance, enabling weight reductions. The value provided in this option is a
        multiplier on the outer wall thickness. The default value of 0.8 is higher than it
        would be if it were purely empirically determined from Sullivan et al. 2006
        (https://ntrs.nasa.gov/citations/20060021606), but has been made much more
        conservative to fall more in line with ~60% gravimetric efficiency tanks
    youngs_modulus : float
        Young's modulus of wall material (Pa), by default LiAl 2090 taken from Table XIII of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    density : float
        Density of wall material (kg/m^3), by default LiAl 2090 taken from Table XIII of
        Sullivan et al. 2006 (https://ntrs.nasa.gov/citations/20060021606)
    """

    def initialize(self):
        self.options.declare("safety_factor", default=2.0, desc="Safety factor on wall thickness")
        self.options.declare("stiffening_multiplier", default=0.8, desc="Multiplier on wall thickness")
        self.options.declare("youngs_modulus", default=8.0e10, desc="Young's modulus of wall material in Pa")
        self.options.declare("density", default=2699.0, desc="Density of wall material in kg/m^3")

    def setup(self):
        self.add_input("design_pressure_differential", val=101325.0, units="Pa")
        self.add_input("radius", val=radius_constraint, units="m") # by default 0.5m, 1m used for the analysis.
        self.add_input("length", val=tank_length, units="m") # by default 2m, 10m used for the analysis. 

        self.add_output("thickness", lower=0.0, units="m")
        self.add_output("weight", lower=0.0, units="kg")

        self.declare_partials(["thickness", "weight"], ["design_pressure_differential", "radius", "length"])

    def compute(self, inputs, outputs):
        p = inputs["design_pressure_differential"]
        r = inputs["radius"]
        L = inputs["length"]
        SF = self.options["safety_factor"]
        E = self.options["youngs_modulus"]
        density = self.options["density"]
        stiff_mult = self.options["stiffening_multiplier"]

        # Compute the thickness necessary for the cylindrical portion
        t_cyl = (p * SF * L * r**1.5 / (0.92 * E)) ** (1 / 2.5)

        # Compute the thickness necessary for the spherical portion
        t_sph = r * np.sqrt(p * SF / (0.365 * E))

        # Take the maximum of the two, when r and L are small the KS
        # isn't a great approximation and the weighting parameter needs
        # to be very high, so just let it be C1 discontinuous
        outputs["thickness"] = stiff_mult * np.maximum(t_cyl, t_sph)

        surface_area = 4 * np.pi * r**2 + 2 * np.pi * r * L
        outputs["weight"] = surface_area * outputs["thickness"] * density

    def compute_partials(self, inputs, J):
        p = inputs["design_pressure_differential"]
        r = inputs["radius"]
        L = inputs["length"]
        SF = self.options["safety_factor"]
        E = self.options["youngs_modulus"]
        density = self.options["density"]
        stiff_mult = self.options["stiffening_multiplier"]

        # Compute the thickness necessary for the cylindrical portion
        t_cyl = (p * SF * L * r**1.5 / (0.92 * E)) ** (1 / 2.5)
        if L < 1e-6:
            dtcyl_dp = 0.0
            dtcyl_dr = 0.0
            dtcyl_dL = 0.0
        else:
            first_term = (p * SF * L * r**1.5 / (0.92 * E)) ** (1 / 2.5 - 1) / 2.5
            dtcyl_dp = first_term * SF * L * r**1.5 / (0.92 * E)
            dtcyl_dr = first_term * p * SF * L * r**0.5 / (0.92 * E) * 1.5
            dtcyl_dL = first_term * p * SF * r**1.5 / (0.92 * E)

        # Compute the thickness necessary for the spherical portion
        t_sph = r * np.sqrt(p * SF / (0.365 * E))
        dtsph_dp = 0.5 * r * (p * SF / (0.365 * E)) ** (-0.5) * SF / (0.365 * E)
        dtsph_dr = t_sph / r
        dtsph_dL = 0.0

        # Derivative is from whichever thickness is greater
        use_cyl = t_cyl.item() > t_sph.item()
        J["thickness", "design_pressure_differential"] = (dtcyl_dp if use_cyl else dtsph_dp) * stiff_mult
        J["thickness", "radius"] = (dtcyl_dr if use_cyl else dtsph_dr) * stiff_mult
        J["thickness", "length"] = (dtcyl_dL if use_cyl else dtsph_dL) * stiff_mult

        t = stiff_mult * np.maximum(t_cyl, t_sph)
        A = 4 * np.pi * r**2 + 2 * np.pi * r * L
        dAdr = 8 * np.pi * r + 2 * np.pi * L
        dAdL = 2 * np.pi * r
        J["weight", "design_pressure_differential"] = A * J["thickness", "design_pressure_differential"] * density
        J["weight", "radius"] = (dAdr * t + A * J["thickness", "radius"]) * density
        J["weight", "length"] = (dAdL * t + A * J["thickness", "length"]) * density

class MLIWeight(om.ExplicitComponent):
    """
    Compute the weight of the MLI given the tank geometry and number of MLI layers.
    Foil and spacer areal density per layer estimated from here:
    https://frakoterm.com/cryogenics/multi-layer-insulation-mli/

    Inputs
    ------
    radius : float
        Inner radius of the cylinder and hemispherical end caps. This value
        does not include the insulation (scalar, m).
    length : float
        Length of JUST THE CYLIDRICAL part of the tank (scalar, m)
    N_layers : float
        Number of reflective sheild layers in the MLI, should be at least ~10 for model
        to retain reasonable accuracy (scalar, dimensionless)

    Outputs
    -------
    weight : float
        Total weight of the MLI insulation (scalar, kg)

    Options
    -------
    foil_layer_areal_weight : float
        Areal weight of a single foil layer, by default 18e-3 (scalar, kg/m^2)
    spacer_layer_areal_weight : float
        Areal weight of a single spacer layer, by default 12e-3 (scalar, kg/m^2)
    """

    def initialize(self):
        self.options.declare("foil_layer_areal_weight", default=18e-3, desc="Areal weight of foil layer in kg/m^2")
        self.options.declare("spacer_layer_areal_weight", default=12e-3, desc="Areal weight of spacer layer in kg/m^2")

    def setup(self):
        self.add_input("radius", units="m")
        self.add_input("length", units="m")
        self.add_input("N_layers")

        self.add_output("weight", units="kg")

        self.declare_partials("weight", ["radius", "length", "N_layers"])

    def compute(self, inputs, outputs):
        r = inputs["radius"]
        L = inputs["length"]
        N = inputs["N_layers"]
        W_foil = self.options["foil_layer_areal_weight"]
        W_spacer = self.options["spacer_layer_areal_weight"]

        # Compute surface area
        A = 4 * np.pi * r**2 + 2 * np.pi * r * L

        outputs["weight"] = (W_foil + W_spacer) * N * A

    def compute_partials(self, inputs, J):
        r = inputs["radius"]
        L = inputs["length"]
        N = inputs["N_layers"]
        W_foil = self.options["foil_layer_areal_weight"]
        W_spacer = self.options["spacer_layer_areal_weight"]

        # Compute surface area
        A = 4 * np.pi * r**2 + 2 * np.pi * r * L

        J["weight", "N_layers"] = (W_foil + W_spacer) * A
        J["weight", "radius"] = (W_foil + W_spacer) * N * (8 * np.pi * r + 2 * np.pi * L)
        J["weight", "length"] = (W_foil + W_spacer) * N * (2 * np.pi * r)

if __name__ == "__main__":
    p = om.Problem()
    p.model.add_subsystem("model", VacuumTankWeight(), promotes=["*"])
    p.setup(force_alloc_complex=True)

    p.set_val("environment_design_pressure", 1.0, units="atm")
    p.set_val("max_expected_operating_pressure", 2, units="bar") # default is 2.5bar but this looks a bit high given the requirement for storing LH2 is not that high. 
    p.set_val("vacuum_gap", MLI_thickness_inch, units="inch") # default 4 inches. Needs to be at least 2 to maintain the accuracy of the results. 
    p.set_val("radius", radius_constraint_ft, units="ft") # default 8.5 / 2. 1m is 3.28084ft. This is the only radius value that affects the final mass.
    p.set_val("length", tank_length_ft, units="ft") #default 0. 5m is 16.4042ft. This is the only length value that affects the final mass. 

    p.run_model()

    # p.check_partials(method="cs", compact_print=True)

    p.model.list_outputs(units=True)

    r = p.get_val("radius", units="m").item()
    L = p.get_val("length", units="m").item()
    W_LH2 = (4 / 3 * np.pi * r**3 + np.pi * r**2 * L) * 70 * 0.95
    W_tank = p.get_val("weight", units="kg").item()
    print(f"\n-------- Approximate gravimetric efficiency: {W_LH2 / (W_LH2 + W_tank) * 100:.1f}% --------")
    
    print("Tank quantity chosen:", tank_quantity)
    print("Tank weight x tank quantity:", W_tank*tank_quantity, "kg")




##########################STABILITY & CONTROL CODE ##########################

#Static and Dynamic Stability & Control Calculations for a Blended Wing Body Aircraft Design 
#Jude Doherty 28/03/2025 last 
#Import Stability and Control Python programme containing functions
from StabilityControlFunctions.py import WeightBalance, StaticStability, DynamicStabilityControl

#Run Weight and Balance Function
WeightBalance()
#Run Static Stability Analysis
StaticStability()
#Run Dynamic Stability and Control Analysis
DynamicStability()
