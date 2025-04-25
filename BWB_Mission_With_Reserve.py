def missionplot(): 
    import numpy as np

    import openmdao.api as om
    from openconcept.utilities import AddSubtractComp, DictIndepVarComp, plot_trajectory

    # imports for the airplane model itself
    from openconcept.aerodynamics import PolarDrag
    from openconcept.examples.aircraft_data.BWB import data as acdata
    from openconcept.mission import MissionWithReserve, IntegratorGroup
    from openconcept.propulsion import RubberizedTurbofan
    from openconcept.utilities import Integrator, AddSubtractComp, ElementMultiplyDivideComp, DictIndepVarComp
    hydrogen = True

    class B738AirplaneModel(IntegratorGroup):
        """
        A custom model specific to the Boeing 737-800 airplane.
        This class will be passed in to the mission analysis code.

        """

        def initialize(self):
            self.options.declare("num_nodes", default=1)
            self.options.declare("flight_phase", default=None)

        def setup(self):
            nn = self.options["num_nodes"]
            flight_phase = self.options["flight_phase"]

            # a propulsion system needs to be defined in order to provide thrust
            # information for the mission analysis code
            # ==============================================================================
            # Propulsion
            # ==============================================================================
            # -------------- CFM56 engine surrogate model --------------
            self.add_subsystem(
                "N3",
                RubberizedTurbofan(num_nodes=nn, hydrogen = True if hydrogen else False, engine = "N3"),
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

            # use a different drag coefficient for takeoff versus cruise
            if flight_phase not in ["v0v1", "v1v0", "v1vr", "rotate"]:
                cd0_source = "ac|aero|polar|CD0_cruise"
            else:
                cd0_source = "ac|aero|polar|CD0_TO"
            self.add_subsystem(
                "drag",
                PolarDrag(num_nodes=nn),
                promotes_inputs=["fltcond|CL", "ac|geom|*", ("CD0", cd0_source), "fltcond|q", ("e", "ac|aero|polar|e")],
                promotes_outputs=["drag"],
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

            # generally the weights module will be custom to each airplane
            passthru = om.ExecComp("OEW=x", x={"val": 1.0, "units": "kg"}, OEW={"val": 1.0, "units": "kg"})
            self.add_subsystem("OEW", passthru, promotes_inputs=[("x", "ac|weights|OEW")], promotes_outputs=["OEW"])

            self.add_subsystem(
                "weight",
                AddSubtractComp(
                    output_name="weight",
                    input_names=["ac|weights|MTOW", "fuel_burn"],
                    units="kg",
                    vec_size=[1, nn],
                    scaling_factors=[1, -1],
                ),
                promotes_inputs=["*"],
                promotes_outputs=["weight"],
            )


    class BWBAnalysisGroup(om.Group):
        def setup(self):
            # Define number of analysis points to run pers mission segment
            nn = 11

            # Define a bunch of design varaiables and airplane-specific parameters
            dv_comp = self.add_subsystem("dv_comp", DictIndepVarComp(acdata), promotes_outputs=["*"])
            dv_comp.add_output_from_dict("ac|aero|CLmax_TO")
            dv_comp.add_output_from_dict("ac|aero|polar|e")
            dv_comp.add_output_from_dict("ac|aero|polar|CD0_TO")
            dv_comp.add_output_from_dict("ac|aero|polar|CD0_cruise")

            dv_comp.add_output_from_dict("ac|geom|wing|S_ref")
            dv_comp.add_output_from_dict("ac|geom|wing|AR")
            dv_comp.add_output_from_dict("ac|geom|wing|c4sweep")
            dv_comp.add_output_from_dict("ac|geom|wing|taper")
            dv_comp.add_output_from_dict("ac|geom|wing|toverc")
            dv_comp.add_output_from_dict("ac|geom|hstab|S_ref")
            dv_comp.add_output_from_dict("ac|geom|hstab|c4_to_wing_c4")
            dv_comp.add_output_from_dict("ac|geom|vstab|S_ref")

            dv_comp.add_output_from_dict("ac|geom|nosegear|length")
            dv_comp.add_output_from_dict("ac|geom|maingear|length")

            dv_comp.add_output_from_dict("ac|weights|MTOW")
            dv_comp.add_output_from_dict("ac|weights|W_fuel_max")
            dv_comp.add_output_from_dict("ac|weights|MLW")
            dv_comp.add_output_from_dict("ac|weights|OEW")

            dv_comp.add_output_from_dict("ac|propulsion|engine|rating")

            dv_comp.add_output_from_dict("ac|num_passengers_max")
            dv_comp.add_output_from_dict("ac|q_cruise")

            # Run a full mission analysis including takeoff, reserve_, cruise,reserve_ and descereserve_nt
            self.add_subsystem(
                "analysis",
                MissionWithReserve(num_nodes=nn, aircraft_model=B738AirplaneModel),
                promotes_inputs=["*"],
                promotes_outputs=["*"],
            )


    def configure_problem():
        prob = om.Problem()
        prob.model = BWBAnalysisGroup()
        prob.model.nonlinear_solver = om.NewtonSolver(iprint=2, solve_subsystems=True)
        prob.model.linear_solver = om.DirectSolver()
        prob.model.nonlinear_solver.options["maxiter"] = 20
        prob.model.nonlinear_solver.options["atol"] = 1e-6
        prob.model.nonlinear_solver.options["rtol"] = 1e-6
        prob.model.nonlinear_solver.linesearch = om.BoundsEnforceLS(bound_enforcement="scalar", print_bound_enforce=False)
        return prob


    def set_values(prob, num_nodes):
        # set some (required) mission parameters. Each pahse needs a vertical and air-speed
        # the entire mission needs a cruise altitude and range
        prob.set_val("climb.fltcond|vs", np.linspace(2300.0, 600.0, num_nodes), units="ft/min")
        prob.set_val("climb.fltcond|Ueas", np.linspace(230, 220, num_nodes), units="kn")
        prob.set_val("cruise.fltcond|vs", np.ones((num_nodes,)) * 0.0, units="ft/min")
        prob.set_val("cruise.fltcond|Ueas", np.linspace(260, 260, num_nodes), units="kn")
        prob.set_val("descent.fltcond|vs", np.linspace(-1000, -150, num_nodes), units="ft/min")
        prob.set_val("descent.fltcond|Ueas", np.ones((num_nodes,)) * 250, units="kn")
        prob.set_val("reserve_climb.fltcond|vs", np.linspace(3000.0, 2300.0, num_nodes), units="ft/min")
        prob.set_val("reserve_climb.fltcond|Ueas", np.linspace(230, 230, num_nodes), units="kn")
        prob.set_val("reserve_cruise.fltcond|vs", np.ones((num_nodes,)) * 0.0, units="ft/min")
        prob.set_val("reserve_cruise.fltcond|Ueas", np.linspace(242, 242, num_nodes), units="kn")
        prob.set_val("reserve_descent.fltcond|vs", np.linspace(-800, -800, num_nodes), units="ft/min")
        prob.set_val("reserve_descent.fltcond|Ueas", np.ones((num_nodes,)) * 250, units="kn")
        prob.set_val("loiter.fltcond|vs", np.linspace(0.0, 0.0, num_nodes), units="ft/min")
        prob.set_val("loiter.fltcond|Ueas", np.ones((num_nodes,)) * 200, units="kn")
        prob.set_val("cruise|h0", 30000.0, units="ft")
        prob.set_val("reserve|h0", 15000.0, units="ft")
        prob.set_val("mission_range", 2000, units="NM")


    def show_outputs(prob):
        # print some outputs
        #vars_list = ["descent.fuel_used_final", "loiter.fuel_used_final"]
        #units = ["lb", "lb"]
        #nice_print_names = ["Block fuel", "Total fuel"]
        #print("=======================================================================")
        #for i, thing in enumerate(vars_list):
        #    print(nice_print_names[i] + ": " + str(prob.get_val(thing, units=units[i])[0]) + " " + units[i])

        # plot some stuff
        plots = True
        if plots:
            x_var = "range"
            x_unit = "NM"
            y_vars = ["fltcond|h", "fltcond|M"]
            y_units = ["ft", None]
            x_label = "Range (nmi)"
            y_labels = [
                "Altitude (ft)",
                "Mach number",   
            ]
            phases = ["climb", "cruise", "descent", "reserve_climb", "reserve_cruise", "reserve_descent", "loiter"]
            plot_trajectory(
                prob,
                x_var,
                x_unit,
                y_vars,
                y_units,
                phases,
                x_label=x_label,
                y_labels=y_labels,
                marker="-",
                plot_title="BWB Mission Profile",
            )


    def run_BWB_analysis(plots=False):
        num_nodes = 11
        prob = configure_problem()
        prob.setup(check=True, mode="fwd")
        set_values(prob, num_nodes)
        prob.run_model()
        prob.model.list_outputs()
        if plots:
            show_outputs(prob)
        return prob


    if __name__ == "__main__":
        run_BWB_analysis(plots=True)

def takeoff_and_landing():
    from scipy.integrate import cumulative_trapezoid
    import matplotlib.pyplot as plt

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import cumulative_trapezoid


    ##################################################################
    #                        Takeoff                                 #
    ##################################################################
    def takeoff_distance(W, rho, S, CL_max, CD0, T, mu, g, V_R, V1, V2, k1, k2, phi):
        
        V = np.linspace(0, V_R, 100)  # Speed range from 0 to rotation speed
        q = 0.5 * rho * V**2  # Dynamic pressure
        L = q * S * CL_max  # Lift force
        CD = CD0 + (phi* k1 * (CL_max**4)) + (phi * k2 * (CL_max**2))
        D = q * S * CD  # Drag force
        R = mu * (W - L)  # Rolling resistance
        F_net = T - (D + R)  # Net accelerating force
        a = g*F_net  # Acceleration
        
        # Integrate acceleration to get velocity squared vs. distance
        V_sq = V**2
        S1 = cumulative_trapezoid(W / (2 * a), V_sq, initial=0)  # Ground roll distance
        
        return S1, V

    def takeoff_distance_reduced_thrust(W, rho, S, CL_max, CD0, T, mu, g, V_start, V_R, S_start):
        T_half = 0.5 * T
        V = np.linspace(V_start, V_R, num=80)
        q = 0.5 * rho * V**2
        L = q * S * CL_max
        D = q * S * CD0
        R = mu * (W - L)
        F_net = T_half - (D + R)
        a = g*F_net
        V_sq = V**2
        S2 = cumulative_trapezoid(W / (2 * a), V_sq, initial=0)
        S2 += S_start  # Start at point where V1 occurs
        return S2, V

    def abort_landing_distance(W, mu_braking, g, V1, S1_array, V):

        a_braking = -mu_braking * g  # Constant deceleration
        S1_final = np.interp(V1, V, S1_array) 
        S_abort = (V1**2) / (2 * abs(a_braking))  # Distance using suvat
        
        V_abort = np.linspace(V1, 0, 50)
        S_abort_array = np.linspace(S1_final, S1_final + S_abort, 50)
        
        return S_abort, V_abort, S_abort_array

    def climb_angle(T, D, W):

        climb_ratio = (T - D) / W
        theta_climb = np.degrees(np.arcsin(climb_ratio))
        return theta_climb


    def obstacle_clearance_distance(climb_angle_deg, V_R, g):
        Rotate = (V_R**2)/(0.5*g)
        S3a = Rotate*np.sin(climb_angle_deg)
        return S3a, Rotate


    def clear35ft(S3a, Rotate, climb_angle_deg):
        h = (S3a**2)/(2*Rotate)
        S3b = (10.7 - h)/0.2788 ###FIX - TAN FUNC NOT WORKING
        return S3b, h

    def ground_effect(h, b):
        numerator = (((16*h)/b)**2)
        denomenator = (1+(((16*h)/b)**2))
        phi = numerator/denomenator 
        return phi

    # parameters
    W = 900000  
    rho = 1.225  
    S = 70 #just based off the wings  
    CL_max = 1.3  
    CD0 = 0.0314  
    T = 248210  
    mu = 0.02  
    g = 9.81     
    mu_braking = 0.4  
    h_obstacle = 35 * 0.3048  #metres
    h = 2 # wing height off ground
    b = 37.638 #wingspan
    k1 = 0.1373
    k2 = 0.0205

    V_Stall = 0.5*np.sqrt((2 * W) / (rho * S * CL_max))
    V_R = 1.2*V_Stall
    V_mc = 0.95*V_Stall
    V1 = V_mc
    V2 = V_R
    V_final = V_R
    phi = ground_effect(h, b)
    S1_array, V = takeoff_distance(W, rho, S, CL_max, CD0, T, mu, g, V_R, V1, V2, k1, k2, phi)
    S1_V1 = np.interp(V1, V, S1_array)
    S_abort, V_abort, S_abort_array = abort_landing_distance(W, mu_braking, g, V1, S1_array, V)
    D = 0.5 * rho * V_R**2 * S * CD0  #drag at rotation speed
    climb_angle_deg = climb_angle(T, D, W)
    S3a, Rotate = obstacle_clearance_distance(climb_angle_deg, V_R, g)
    S3b, h = clear35ft(S3a, Rotate, climb_angle_deg)
    clear_obstacle = S3a + S3b
    Stot = S1_array[-1] + S3a + S3b  # Total distance including ground roll
    S2_array, V_reduced = takeoff_distance_reduced_thrust(W, rho, S, CL_max, CD0, T, mu, g, V1, V_final, S1_V1)
    total_ground_roll = S1_array[-1]

    print("Takeoff")
    print(f"Estimated Ground Roll Distance at V1: {np.interp(V1, V, S1_array):.2f} meters")
    print(f"Estimated Aborted Takeoff Distance: {S_abort:.2f} meters")
    print(f"Estimated Climb Angle: {climb_angle_deg:.2f} degrees")
    print(f"Estimated Stall Speed: {V_Stall:.2f} m/s")
    print(f"Estimated Distance to Clear 35ft Obstacle (S3a + S3b): {clear_obstacle:.2f} meters")
    print(f"Total Takeoff Distance: {Stot:.2f} meters")
    print(f"Total Ground Roll Distance: {total_ground_roll:.2f} meters")

    ##################################################################
    #                        Landing                                 #
    ##################################################################

    def CL_Steady_Approach(CL_max):
        """Calculate the airborne distance SA."""
        CL_A = CL_max/1.69
        return CL_A

    def airborne_distance(theta_D, g, V_A, R_flare):
        """Calculate the airborne distance SA."""
        SA = (15.2 / (theta_D)) + ((R_flare)*(theta_D) / (2))
        return SA

    def flare_radius(V_A, CL_A, CL_R, g):
        """Calculate the flare radius R."""
        R_flare = (V_A ** 2) / (g * (CL_R / CL_A - 1))
        return R_flare

    def transition_distance(V_A):
        """Calculate the transition distance."""
        S_tran = 2 * V_A
        return S_tran

    def calculate_landing_ground_roll(W_landing, S, C_D0, k1, k2, CL_max, mu, V_A, rho, g, CL_A):    
        # array of velocities from 0 to V_A
        V = np.linspace(0.00001, V_A, 100)
        V_squared = V**2
        # Calculate CL during ground roll (assume constant angle of attack)
        CL = CL_A * (V_A**2 / V_squared)
        CL = np.where(V > 0.1, CL, CL_max)
        CD = C_D0 + (k2 * CL**2) + (k1 * CL**4)
        L = 0.5 * rho * V**2 * S * CL
        D = 0.5 * rho * V**2 * S * CD
        # Calculate the integrand: W / (2g * [T_R + D + μ(W - L)])
        denominator = ((D + mu * (W - L)))
        integrand = W_landing / (2 * g * denominator)
        # Handle division by zero at V=0 (set to integrand at first non-zero velocity)
        integrand[0] = integrand[1]
        ground_roll_distance = np.trapz(integrand, V_squared)
        
        return ground_roll_distance, V_squared, integrand


    # Aircraft parameters
    W_landing = 500000            
    C_D0 = 0.0314      
    k1 = 0.1373
    k2 = 0.0205            
    V_A = 70       
    theta_D = 3  
    CL_A = CL_Steady_Approach(CL_max)  
    CL_R = 2  
        
    # ground roll distance
    distance, V_sq, integrand = calculate_landing_ground_roll(W_landing, S, C_D0, k1, k2, CL_max, mu, V_A, rho, g, CL_A)

    R_flare = flare_radius(V_A, CL_A, CL_R, g)
    SA = airborne_distance(theta_D, g, V_A, R_flare)
    S_tran = transition_distance(V_A)

    total_distance = [distance, SA, S_tran]
    landing = sum(total_distance)
    print("Landing")
    print(f"Approach Speed: {V_A:.1f} m/s ({V_A * 1.944:.1f} knots)")
    print(f"Airbourne distance from altitude of 15.2m to touchdown: {SA:.1f} metres")
    print(f"Total Flare radius: {R_flare:.1f} metres")
    print(f"Total transition distance from landing to braking: {S_tran:.1f} metres")
    print(f"Landing Ground Roll Distance: {distance:.1f} meters")
    print(f"Total approach and landing distance: {landing:.1f} metres")

    ################################################
    #                   Plot                       #
    ################################################
    plt.plot(S1_array, V, label="Ground Roll Velocity")
    plt.plot(S2_array, V_reduced, 'b--', label="Reduced Thrust (50%) Segment")
    plt.axhline(V1, color='r', linestyle='--', label="V1 (Decision Speed)")
    plt.axhline(V2, color='g', linestyle='--', label="V2 (Takeoff Safety Speed)")
    plt.plot(S_abort_array, V_abort, 'm--', label="Abort Deceleration")
    plt.ylabel("Velocity (m/s)")
    plt.xlabel("Distance (m)")
    plt.title("Aircraft Takeoff and Aborted Takeoff Velocities")
    plt.legend()
    plt.grid()
    plt.show()

missionplot()
takeoff_and_landing()