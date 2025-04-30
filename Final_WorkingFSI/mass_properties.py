def write_aircraft_data_to_csv(filename, cabin_mass_array, wing_mass_array, vstab_mass_array, landing_gear_mass_array):

    fuselage_mass = cabin_mass_array[0]
    fuselage_x_cg = cabin_mass_array[1]
    fuselage_y_cg = cabin_mass_array[2]
    fuselage_z_cg = cabin_mass_array[3]
    fuselage_Ixx = cabin_mass_array[4]
    fuselage_Iyy = cabin_mass_array[5]
    fuselage_Izz = cabin_mass_array[6]

    aft_body_mass = 2*1406
    aft_body_Ixx = 2*2.98836e4
    aft_body_Iyy = 2*8.53539e+03
    aft_body_Izz = 2*1.50423e+04
    aft_body_x_cg = -12.4-1.776
    aft_body_y_cg = 0
    aft_body_z_cg = 2.3-1.7

    wing_x_cg = wing_mass_array[0]
    wing_y_cg = wing_mass_array[1]
    wing_z_cg = wing_mass_array[2]
    wing_Ixx = wing_mass_array[3]
    wing_Iyy = wing_mass_array[4]
    wing_Izz = wing_mass_array[5]
    wing_mass = wing_mass_array[6]
    buckled = wing_mass_array[7]
    if buckled == 1:
        print("WING BUCKLED, ADJUST WING PARAMETERS")

    vstab_x_cg = vstab_mass_array[0]
    vstab_y_cg = vstab_mass_array[1]
    vstab_z_cg = vstab_mass_array[2]
    vstab_Ixx = vstab_mass_array[3]
    vstab_Iyy = vstab_mass_array[4]
    vstab_Izz = vstab_mass_array[5]
    vstab_mass = vstab_mass_array[6]
    buckled = vstab_mass_array[7]
    if buckled == 1:
        print("VERTICAL STABILISER BUCKLED, ADJUST PARAMETERS")

    main_lg_mass = landing_gear_mass_array[0]
    nose_lg_mass = landing_gear_mass_array[1]
    main_lg_cg_from_hinge = landing_gear_mass_array[2]
    nose_lg_cg_from_hinge = landing_gear_mass_array[3]
    main_lg_Ixx = 0
    main_lg_Iyy = 0
    main_lg_Izz = 0
    main_lg_x_cg = 0
    main_lg_y_cg = 0
    main_lg_z_cg = 0

    nose_lg_Ixx = 0
    nose_lg_Iyy = 0
    nose_lg_Izz = 0
    nose_lg_x_cg = 0
    nose_lg_y_cg = 0
    nose_lg_z_cg = 0

    # Data to be written
    data = [
        ['Component', 'mass', 'Ixx', 'Iyy', 'Izz', 'x_cg', 'y_cg', 'z_cg', 'hinge_cg'],
        ['Fuselage', fuselage_mass, fuselage_Ixx, fuselage_Iyy, fuselage_Izz, fuselage_x_cg, fuselage_y_cg, fuselage_z_cg, 0],
        ['Aft Body', aft_body_mass, aft_body_Ixx, aft_body_Iyy, aft_body_Izz, aft_body_x_cg, aft_body_y_cg, aft_body_z_cg, 0],
        ['Wing', wing_mass, wing_Ixx, wing_Iyy, wing_Izz, wing_x_cg, wing_y_cg, wing_z_cg, 0],
        ['Main Landing Gear', main_lg_mass, main_lg_Ixx, main_lg_Iyy, main_lg_Izz, main_lg_x_cg, main_lg_y_cg, main_lg_z_cg, main_lg_cg_from_hinge],
        ['Nose Landing Gear', nose_lg_mass, nose_lg_Ixx, nose_lg_Iyy, nose_lg_Izz, nose_lg_x_cg, nose_lg_y_cg, nose_lg_z_cg, nose_lg_cg_from_hinge],
        ['vstab', vstab_mass, vstab_Ixx, vstab_Iyy, vstab_Izz, vstab_x_cg, vstab_y_cg, vstab_z_cg, 0],
    ]
    import csv
    # Write to CSV
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)