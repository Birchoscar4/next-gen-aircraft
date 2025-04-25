import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt

def calculate_bending_moment(loads, rib_pitch):
    n = len(loads)
    bending_moment = np.zeros(n)

    # Compute bending moment at each force application point
    for i in range(n):
        moment = 0
        for j in range(i, n):
            moment += loads[j] * ((j - i) * rib_pitch)
        bending_moment[i] = moment

    return bending_moment

xpos = np.linspace(0, 10, 11)
loads = np.linspace(1000, 1000, 11)

bending_moment=calculate_bending_moment(loads, rib_pitch=1.2)

#plt.plot(xpos, bending_moment)
#plt.show()

def wingbox_polar_moment(spar_thickness, wingbox_height, skin_thickness, wingbox_chord, I_spar_cap, I_stringer, A_spar_cap, A_stringer, num_stringers):
    # Contribution of the spars to Ix (bending about z-axis)
    I_spar_x = (spar_thickness * (wingbox_height ** 3)) / 12
    
    # Contribution of the spars to Iy (bending about y-axis)
    I_spar_y = (wingbox_height * (spar_thickness ** 3)) / 12 + (wingbox_height * spar_thickness) * (wingbox_chord / 2) ** 2
    
    # Contribution of the top and bottom skins to Ix
    I_skin_x = (wingbox_chord * (skin_thickness ** 3)) / 12 + (wingbox_chord * skin_thickness) * (wingbox_height / 2) ** 2
    
    # Contribution of the top and bottom skins to Iy
    I_skin_y = (skin_thickness * (wingbox_chord ** 3)) / 12 
    
    # Contribution of stringers (if present)
    if num_stringers > 0:
        I_stringers_x = I_stringer + A_stringer * (wingbox_height / 2) ** 2
        I_stringers_y = I_stringer + A_stringer * (wingbox_chord / 2) ** 2
    else:
        I_stringers_x = 0
        I_stringers_y = 0
    
    # Spar caps at the corners of the wingbox
    I_spar_caps_x = I_spar_cap + A_spar_cap * (wingbox_height / 2) ** 2
    I_spar_caps_y = I_spar_cap + A_spar_cap * (wingbox_chord / 2) ** 2
    
    # Total second moment of area about x and y axes
    I_x_total = 2 * I_spar_x + 4 * I_spar_caps_x + 2 * I_skin_x + 2 * num_stringers * I_stringers_x
    I_y_total = 2 * I_spar_y + 4 * I_spar_caps_y + 2 * I_skin_y + 2 * num_stringers * I_stringers_y
    
    # Polar moment of inertia (J = Ix + Iy)
    J_total = I_x_total + I_y_total
    
    return J_total

polar_moment = wingbox_polar_moment(0.01,0.28,0.002,1,0,0,0,0,0)
print ("Polar moment",polar_moment)

def skin_second_moment(skin_thickness, wingbox_chord, I_stringer, A_stringer, num_stringers):
    I_skin = (wingbox_chord * (skin_thickness ** 3)) / 12
    # Contribution of stringers
    if num_stringers > 0:
        I_stringer =  I_stringer 
    else:
        I_stringer = 0

    # Total second moment of area
    I_skin_total =  I_skin + num_stringers * I_stringer
    return I_skin_total

#skinEI = 70e9 * skin_second_moment(0.02, 1, )

def wingbox_polar_moment(spar_thickness, wingbox_height, skin_thickness, wingbox_chord):

    A = wingbox_height * wingbox_chord

    # Lengths of each section of the closed wingbox
    a = wingbox_chord  # top and bottom skins
    b = wingbox_height  # front and rear spars
    t1 = skin_thickness
    t2 = spar_thickness

    J = (2*t1*t2*(a-t2)**2*(b-t1)**2)/(a*t2+b*t1-t1**2-t2**2)
    print(a*t1)
    print(b*t2)
    print(t1**2)
    print(t2**2)
    print(a*t1+b*t2-t1**2-t2**2)
    return J

print(wingbox_polar_moment(4,200,2,1000))