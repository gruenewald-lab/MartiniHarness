"""
Minimal polyply/GROMACS simulation setup and execution script.
"""

import logging
from pathlib import Path
import time
import os
import gromacs
from polyply import gen_coords
from polyply.src.topology import Topology
from vermouth.forcefield import ForceField
import numpy as np

log = logging.getLogger(__name__)

base_dir=Path("/home/fgruenewald/ProgramDev/MartiniHarness")
ff_dir=base_dir/"martini_harness/data/gmx_files"
mdp_dir=base_dir/"martini_harness/data/mdp_files")
FORCEFIELD_ITP   = ff_dir / Path("martini_v3.0.0.itp")
SOLVENT_ITP   = ff_dir / Path("martini_v3.0.0_solvents_v1.itp")
MOLECULE_ITP_DIR = base_dir / "martini_harness/data/row42.itp"
MIN_MDP          = mdp_dir / Path("min.mdp")
EQ_MDP           = mdp_dir  / Path("eq.mdp")


def select_solvent_and_target(molecules: list[str]) -> tuple[str, str]:
    solvent, target = np.random.choice(molecules, size=2)
    return solvent, target

def write_topology(solvent, target, n_target, n_solvent, out_top):
    lines = [f'#include "{FORCEFIELD_ITP}"', f'#include "{SOLVENT_ITP}"']
    lines.append(f'#include "{MOLECULE_ITP_DIR}"')
    lines += [
        "", "[ system ]", f"{target}_in_{solvent}",
        "", "[ molecules ]",
        f"{target}  {n_target}",
        f"{solvent}  {n_solvent}",
    ]
    out_top.write_text("\n".join(lines))

def run_simulation(
    target: str,
    solvent: str,
    n_target: int = 1,
    n_solvent: int = 8000,
    density: float = 900.0,
    pforce: float = 800000.0,
    workdir: Path = Path("simulation"),
):
    print("Making workspace")
    workdir.mkdir(parents=True, exist_ok=True)

    top = workdir / "system.top"
    write_topology(solvent, target, n_target, n_solvent, top)
    print("Generating coordiantes")
    # Generate coordinates with polyply
    gen_coords(toppath=str(top), outpath= workdir / "init.gro", density=density, name="coords")


    print("Running Minimization")
    # Minimisation
    gromacs.grompp(f=MIN_MDP, c=workdir / "init.gro", p=top, o=workdir / "min.tpr", maxwarn=4)
    try:
        gromacs.mdrun(s=workdir / "min.tpr", deffnm=workdir / "min", v=True, nsteps=5000)
    except Exception as e:
        log.warning("Minimisation mdrun failed (possible numerical issue): %s", e)

    # Equilibration
    min_gro = workdir / "min.gro"
    if not min_gro.exists():
        log.warning("min.gro not found, falling back to init.gro")
        min_gro = workdir / "init.gro"

    gromacs.grompp(f=EQ_MDP, c=min_gro, p=top, o=workdir / "eq.tpr", maxwarn = 4)
    try:
        gromacs.mdrun(s=workdir / "eq.tpr", deffnm=workdir / "eq", pforce=pforce, v=True)
    except Exception as e:
        log.warning("Equilibration mdrun failed (possible numerical issue): %s", e)

    # prodiction
    eq_gro = workdir / "eq.gro"
    if not eq_gro.exists():
        raise IOError("eq.gro not found")

    with open(mdp_dir/"awh.mdp") as f:
        text = f.read()

    with open(workdir / "awh.mdp", "w") as f:
        f.write(text.replace("TOCOUPLE", target))

    gromacs.grompp(f=workdir/"awh.mdp", c=eq_gro, p=top, o=workdir / "awh.tpr", maxwarn = 4)
    try:
        gromacs.mdrun(s=workdir / "awh.tpr", deffnm=workdir / "awh", pforce=pforce, v=True)
    except Exception as e:
        log.warning("Equilibration mdrun failed (possible numerical issue): %s", e)

ff = ForceField("test")
top=Topology.from_gmx_topfile("row42.itp", name='test')
mol_names = list(top.force_field.blocks.keys())
for target in mol_names:
    for solvent in ["OCO", "W", "CLF"]:
        cwd = os.getcwd()
        tag = f"{target}_in_{solvent}"
        run_simulation(target, solvent, workdir=cwd / Path(tag))
