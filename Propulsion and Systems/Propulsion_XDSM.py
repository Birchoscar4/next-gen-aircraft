
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


x.add_system("solver", SOLVER, (r"\text{0, 3} \rightarrow \text{1:}", r"\text{Propulsion}"))
x.add_system("enginecycle", FUNC, (r"\text{1:}", r"\text{Engine Cycle}"))
x.add_system("lh2tanksize", FUNC, (r"\text{2:}", r"\text{LH}_2 \text{Tank Size}" ))
x.add_system("lh2tankweight", FUNC, (r"\text{3:}", r"\text{LH}_2 \text{Tank Weight}" ))

x.add_input("enginecycle", ("MTOW_{estimate}", r"\text{BWB aerodynamic parameters}"))
x.add_input("lh2tanksize", (r"\text{Tank shape}", r"\text{Tank configuration}"))
x.add_input("lh2tankweight", (r"\text{Tank material}", r"\text{Insulation method}"))

x.add_process(
    ["solver", "enginecycle", "lh2tanksize", "lh2tankweight"],
    arrow=True
)

x.connect("solver", "enginecycle", ("mdot_{fuel}", "T_{required}"))
x.connect("enginecycle", "lh2tanksize", (r"\text{m}_{LH_2}"))
x.connect("lh2tanksize", "lh2tankweight", ("L_{tank}", "r_{tank}"))

x.connect("enginecycle", "solver", ("BPR", r"\text{LPT and HPT efficiencies}", r"\text{LPC and HPC efficiencies}", "FPR"))
x.connect("lh2tanksize", "solver", (r"\text{V}_{LH_2}", r"\text{Tank quantity}", "L_{available}", "w_{available}"))
x.connect("lh2tankweight", "solver", (r"\text{Tank quantity}", "t_{insulation}"))

x.add_output("enginecycle", ("mdot_{fuel}^*", "T_{required}^*", "MTOW_{updated}"), side=LEFT);
x.add_output("lh2tanksize",("L_{tank}^*", "r_{tank}^*"), side=LEFT);
x.add_output("lh2tankweight", ("m_{tank}^*", "GE_{tank}"), side=LEFT);

x.write("XDSM-Propulsion-final-v2")
