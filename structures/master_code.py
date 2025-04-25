#Aircraft Sizing Loop Master Code Draft 1 

#Last update 02/04/2025  

#################LIST REQUIRED MODULES AND APIS HERE####################### 

import pandas as pd 

import math 

from scipy.signal import TransferFunction, step, lsim 

import numpy as np 

import matplotlib.pyplot as plt 

pi = math.pi 

plt.ion() 

 

 

 

 

 

#####################GEOMETRY AND CONFIGURATION CODE ##################### 

 

 

 

 

 

##########################AERODYNAMICS CODE ############################# 

 

 

 

 

 

 

##########################MISSION ANALYSIS CODE ############################# 

 

 

 

 

 

 

##########################STRUCTURES CODE ############################# 
import csv
from cabin_mass import compute_cabin_mass
from wingbuckling import compute_wing_mass
from landing_gear_sizing import landing_gear_sizing
from mass_properties import write_aircraft_data_to_csv

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
cabin_mass_array = compute_cabin_mass()
# aft_body_mass_array = compute_aft_body_mass()

# Wing bending analysis – FSI HERE 
wing_mass_array = compute_cabin_mass()

# Landing gear sizing – uses MTOW estimate for initial sizing 
landing_gear_mass_array = landing_gear_sizing()

# Write mass properties csv 
write_aircraft_data_to_csv('aircraft_masses.csv',cabin_mass_array, wing_mass_array, landing_gear_mass_array)
 

 

 

 

 

##########################PROP & SYSTEMS CODE ############################# 

 

 

 

 

 

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

 

 

 

 

 

 

 

 