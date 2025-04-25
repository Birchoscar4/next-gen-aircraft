data = dict()
ac = dict()
# ==AERO==================================
aero = dict()
aero["CLmax_TO"] = {"value": 2.0}

polar = dict()
polar["e"] = {"value": 0.82}
polar["CD0_TO"] = {"value": 0.03}
polar["CD0_cruise"] = {"value": 0.01925}

aero["polar"] = polar
ac["aero"] = aero

# ==GEOMETRY==============================
geom = dict()
wing = dict()
wing["S_ref"] = {"value": 28.22, "units": "m**2"}
wing["AR"] = {"value": 3.77}
wing["c4sweep"] = {"value": 25.0, "units": "deg"}
wing["taper"] = {"value": 0.5102}
wing["toverc"] = {"value": 0.14}
geom["wing"] = wing

fuselage = dict()
fuselage["S_wet"] = {"value": 349.21, "units": "ft**2"}
fuselage["width"] = {"value": 49.2, "units": "ft"}
fuselage["length"] = {"value": 58.26, "units": "ft"}
fuselage["height"] = {"value": 2.1312, "units": "ft"}
geom["fuselage"] = fuselage

hstab = dict()
hstab["S_ref"] = {"value": 32.78, "units": "m**2"}
hstab["c4_to_wing_c4"] = {"value": 17.9, "units": "m"}
geom["hstab"] = hstab

vstab = dict()
vstab["S_ref"] = {"value": 12, "units": "m**2"}
geom["vstab"] = vstab

nosegear = dict()
nosegear["length"] = {"value": 3, "units": "ft"}
geom["nosegear"] = nosegear

maingear = dict()
maingear["length"] = {"value": 4, "units": "ft"}
geom["maingear"] = maingear

ac["geom"] = geom

# ==WEIGHTS========================
weights = dict()
weights["MTOW"] = {"value": 80000, "units": "kg"}
weights["OEW"] = {"value": 0.530 * 80000, "units": "kg"}
weights["W_fuel_max"] = {"value": 3024.55, "units": "kg"}
weights["MLW"] = {"value": 66349, "units": "kg"}

ac["weights"] = weights

# ==PROPULSION=====================
propulsion = dict()
engine = dict()
engine["rating"] = {"value": 27900, "units": "lbf"}
propulsion["engine"] = engine

propeller = dict()
propeller["diameter"] = {"value": 2.28, "units": "m"}
propulsion["propeller"] = propeller

ac["propulsion"] = propulsion

# Some additional parameters needed by the empirical weights tools
ac["num_passengers_max"] = {"value": 130}
ac["q_cruise"] = {"value": 212.662, "units": "lb*ft**-2"}
data["ac"] = ac
ac["num_engines"] = {"value": 2}

