
from pyxdsm.XDSM import (
    XDSM,
    OPT,
    SUBOPT,
    SOLVER,
    DOE,
    IFUNC,
    FUNC,
    GROUP,
    IGROUP,
    METAMODEL,
    LEFT,
    RIGHT,
)

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True, optional_latex_packages="soul")


x.add_system("solver", SOLVER, (r"\text{0, 3} \rightarrow \text{1:}", r"\text{Systems}"))
x.add_system("epds", FUNC, (r"\text{1:}", r"\text{Electrical Power Distribution System}", r"\text{(EPDS)}"))
x.add_system("ehs", FUNC, (r"\text{2:}", r"\text{Electro-Hydraulic System}", r"\text{(EHS)}"))
x.add_system("ecs", FUNC, (r"\text{3:}", r"\text{Environmental Control System}", r"\text{(ECS)}"))

x.add_input("epds", (r"\text{Electrical loads required}"))
x.add_input("ehs", (r"\text{Control surface and component quantity}", r"\text{Control surface and component sizes}"))
x.add_input("ecs", (r"\text{Engine heat loads}", r"\text{Aircraft velocity}", r"\text{Number of passengers}"))

x.add_process(
    ["solver", "epds", "ehs", "ecs"],
    arrow=True
)

x.connect("solver", "epds", (r"\text{Electrical loads required}"))
x.connect("epds", "ehs", (r"\text{Electrical subsystem quantity}"))
x.connect("ehs", "ecs", (r"\text{Electrical subsystem heat loads}"))

x.connect("epds", "solver", ("BPR", r"\text{ICAO redundancy requirements}"))
x.connect("ehs", "solver", (r"\text{ICAO redundancy requirements}"))
x.connect("ecs", "solver", (r"\text{Size of heat exchangers required}"))

x.add_output("epds", (r"\text{Electrical component quantity}", r"\text{Total EPDS weight}"), side=LEFT)
x.add_output("ehs", (r"\text{Actuator quantity}", r"\text{Actuator size and specification}", r"\text{Total EHS weight}"), side=LEFT)
x.add_output("ecs", (r"\text{Total ECS weight}"), side=LEFT)

x.write("XDSM-Systems-final")
