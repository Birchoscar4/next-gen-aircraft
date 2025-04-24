"""
@File    :   B738_sizing.py - Edited by OB to consider the BWB values and the 130 passenger requirement as given by the mission specification.
@Date    :   2023/03/25
@Author  :   Eytan Adler
@Description : Data needed B738_sizing OpenConcept example. The data is from a combination of:
- Technical site: http://www.b737.org.uk/techspecsdetailed.htm
- Wikipedia: https://en.wikipedia.org/wiki/Boeing_737#Specifications
- OpenConcept B738 model
"""


data = {
    "ac": {
        # ==============================================================================
        # Aerodynamics
        # ==============================================================================
        "aero": {
            "polar": {
                "e": {"value": 0.82},  # 0.82, given from the Aero team
            },
            "Mach_max": {"value": 0.85},
            "Vstall_land": {"value": 157, "units": "kn"},  # estimate
            "airfoil_Cl_max": {"value": 1.75},  # estimate for supercritical airfoil
            "takeoff_flap_deg": {"value": 15, "units": "deg"},
        },
        # ==============================================================================
        # Propulsion
        # ==============================================================================
        "propulsion": {
            "engine": {
                "rating": {"value": 27e3, "units": "lbf"},
            },
            "num_engines": {"value": 2},
        },
        # ==============================================================================
        # Geometry
        # ==============================================================================
        "geom": {
            # -------------- Wing --------------
            "wing": {
                "S_ref": {"value": 339.545, "units": "m**2"}, # 339.545 from the aero team, 124.6 default.
                "AR": {"value": 5.7978}, # 5.7978 from the aero team, 9.45 default. 
                "c4sweep": {"value": 23.9, "units": "deg"},  # 23.9 from the kink point of the wing from the aero team, 25 by default. 
                "taper": {"value": 0.51}, #0.51 from the aero team, 0.159 by default. 
                "toverc": {"value": 0.14},  # 0.14 from the aero team, 0.12 by default. 
            },
            # -------------- Horizontal stabilizer --------------
            "hstab": {
                # "S_ref": {"value": 32.78, "units": "m**2"},  # not needed since tail volume coefficients are used
                "AR": {"value": 6.16},
                "c4sweep": {"value": 30, "units": "deg"},
                "taper": {"value": 0.203},
                "toverc": {"value": 0.12},  # guess
            },
            # -------------- Vertical stabilizer --------------
            "vstab": {
                # "S_ref": {"value": 26.44, "units": "m**2"},  # not needed since tail volume coefficients are used
                "AR": {"value": 3.2}, # 3.2 from the aero team, 1.91 by default.
                "c4sweep": {"value": 9.8, "units": "deg"}, # 9.8 from the aero team, 35 by default.
                "taper": {"value": 0.25}, # 0.25 from the aero team, 0.271 by default. 
                "toverc": {"value": 0.1},  # 0.1 from the aero team, 0.12 by default.
            },
            # -------------- Fuselage --------------
            "fuselage": {
                "length": {"value": 25.632, "units": "m"}, # 25.632m from the aero team, 38.08 by default.
                "height": {"value": 3.0784, "units": "m"}, # 25.632 * 0.12 thickness from the aero team (3.0784m), 3.76 by default.
            },
            # -------------- Nacelle --------------
            "nacelle": {
                "length": {"value": 2.286, "units": "m"},  # estimate from B787 CFM56 engine
                "diameter": {"value": 4.3, "units": "m"},  # minimum prop diameter required from the mission specification
            },
            # -------------- Main landing gear --------------
            "maingear": {
                "length": {"value": 1.8, "units": "m"},
                "num_wheels": {"value": 6}, # triple bogey configuration as per mission specification
                "num_shock_struts": {"value": 2},
            },
            # -------------- Nose landing gear --------------
            "nosegear": {
                "length": {"value": 1.3, "units": "m"},
                "num_wheels": {"value": 2},
            },
        },
        # ==============================================================================
        # Weights
        # ==============================================================================
        "weights": {
            "W_payload": {"value": 13e3, "units": "kg"}, #this is used to increase the outputted MTOW of the aircraft. For second iteration MTOW needs to be around 80000kg or 176,369lbs. 
        },
        # ==============================================================================
        # Miscellaneous
        # ==============================================================================
        "num_passengers_max": {"value": 130},
        "num_flight_deck_crew": {"value": 2},
        "num_cabin_crew": {"value": 4},
        "cabin_pressure": {"value": 8.95, "units": "psi"},
    },
}
