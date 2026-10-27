# Martini Harness

Run campaigns of Martini simulations defined in a single YAML file. A build
matrix (as in GitHub Actions workflows) expands into one simulation per
combination of target molecule, solvent and parameters. Each simulation
builds a topology, generates coordinates with
[polyply](https://github.com/marrink-lab/polyply_1.0) and runs a configurable
pipeline of GROMACS stages via
[GromacsWrapper](https://github.com/Becksteinlab/GromacsWrapper).

Simulations run one after another. Progress is recorded per simulation, so
rerunning the same command resumes where it stopped and reruns only what
changed.

## Installation

```bash
pip install -e ".[dev]"
```

GROMACS (`gmx`) must be on your `PATH`. Requires Python 3.10+.

## Usage

```bash
martini_harness run examples/sample_campaign.yaml --dry-run   # show what would run
martini_harness run examples/sample_campaign.yaml             # run / resume
martini_harness status examples/sample_campaign.yaml          # progress per step
martini_harness check examples/sample_campaign.yaml           # how many production runs finished
martini_harness analyze examples/sample_campaign.yaml         # run analyses, tabulate results
martini_harness list-molecules examples/row42.itp             # molecule types in ITPs
martini_harness split-itp big.itp -o library/                  # one ITP per molecule
martini_harness sample library/ -n 500 --seed 1 -o targets.yaml  # random subset
```

`run` and `status` accept `--only NAME` (repeatable) to select simulations and
`-o DIR` to override the output directory. `run --force` ignores saved state
and reruns everything. `python -m martini_harness` works as well.

## Campaign file

```yaml
title: Sample Martini solvation campaign
output_dir: ../campaigns/sample        # relative to this file

topology:
  includes:                            # relative to forcefield_dir
    - martini_v3.0.0.itp
    - martini_v3.0.0_solvents_v1.itp

molecules:                             # targets not in the includes
  row42: row42.itp

simulation:                            # defaults, overridable per combination
  name: "{target}_in_{solvent}"
  n_target: 1
  n_solvent: 8000
  density: 900.0

matrix:
  target: [row42]
  solvent: [W, OCO, CLF]
  include:
    - {solvent: W, density: 997.0}
    - {solvent: CLF, n_solvent: 6000, density: 1480.0}
  exclude: []

stages:
  - name: min
    mdp: min.mdp                       # relative to mdp_dir
  - name: eq
    mdp: eq.mdp
    maxwarn: 2
  - name: awh
    mdp: awh.mdp
    substitute: {TOCOUPLE: "{target}"}
```

See [examples/sample_campaign.yaml](examples/sample_campaign.yaml) for a
complete, runnable file.

### Paths

Relative paths are resolved against the directory of the file that sets
them (the campaign file or an [imported file](#imports)), except
`stages[].mdp` (against `mdp_dir`) and `topology.includes` (against
`forcefield_dir`). `forcefield_dir` and `mdp_dir` default to the Martini 3
ITPs and MDP files shipped in `martini_harness/data`.

### Imports

Parts of a campaign, such as a stage pipeline, a molecule library, topology
includes or simulation defaults, can live in separate YAML files and be
reused across campaigns:

```yaml
# campaign.yaml
imports:
  - shared/molecules.yaml            # relative to this file
  - shared/solvation_pipeline.yaml
matrix:
  target: [row42]
  solvent: [W, OCO]
stage_overrides:                     # tweak one stage of the imported pipeline
  awh:
    mdp_params: {nsteps: 500000}
    mdrun: {nt: 8}

# shared/molecules.yaml
molecules:
  row42: ../row42.itp                # relative to shared/

# shared/solvation_pipeline.yaml
stages:
  - {name: min, mdp: min.mdp}
  - {name: eq, mdp: eq.mdp, maxwarn: 2}
  - {name: awh, mdp: awh.mdp, substitute: {TOCOUPLE: "{target}"}}
```

- An imported file is a partial campaign with the same keys. It may import
  further files itself; only the merged result has to be a complete campaign.
- Imports are merged in the order listed, then the importing file's own keys
  are merged on top, so later files win and the importing file wins over
  everything it imports.
- Mappings (`simulation`, `molecules`, `matrix`, ...) are merged key by key.
  Lists (`stages`, `topology.includes`, matrix axes, `include`/`exclude`) and
  scalars are replaced as a whole.
- `stage_overrides` maps a stage name to keys that are merged into that stage
  (mappings key by key, anything else replaced), after all imports. Naming a
  stage that does not exist is an error. It also works on the default
  pipeline.
- An import cycle (`a.yaml` → `b.yaml` → `a.yaml`) or a missing imported file
  is an error. Importing the same file through two routes is fine.

See [examples/imports_campaign.yaml](examples/imports_campaign.yaml), which
builds the sample campaign from the files in `examples/shared/`.

### Matrix

The matrix follows GitHub Actions semantics:

1. Every key other than `include`/`exclude` is an axis; the Cartesian product
   of the axes gives the original combinations.
2. `exclude` entries remove every combination they match (partial matches
   allowed).
3. `include` entries are merged into every original combination whose axis
   values they do not contradict; values added by an earlier include may be
   overwritten. An entry that matches no combination is added as a new one.
   A matrix may consist of `include` entries only.

Every combination needs `target` and `solvent`. Keys named like a
`simulation` field (`name`, `n_target`, `n_solvent`, `density`) override that
field; all keys are available as template variables, e.g. `"{temperature}"`.
Target and solvent must be defined in the topology includes or in
`molecules`; this is checked before anything runs.

### Stages

Every simulation starts with two fixed steps, `topology` (writes
`system.top`) and `coords` (polyply writes `init.gro`), followed by the
stages in order. Each stage runs `grompp` and `mdrun -deffnm <name>`.

| key          | meaning                                                                 |
|--------------|-------------------------------------------------------------------------|
| `name`       | Stage name, used as file prefix (`topology`/`coords` are reserved)      |
| `mdp`        | MDP template                                                            |
| `input`      | Step whose `.gro` is the input structure (default: previous step)       |
| `maxwarn`    | Passed to grompp (default 0)                                            |
| `substitute` | Literal placeholder replacements in the MDP, e.g. `TOCOUPLE: "{target}"`|
| `mdp_params` | MDP parameters to set or override, e.g. `{nsteps: 1000}`                |
| `grompp`     | Extra grompp flags, e.g. `{n: /abs/path/index.ndx}`                     |
| `mdrun`      | Extra mdrun flags, e.g. `{nt: 8, v: true}`                              |

String values may use template variables. If `stages` is omitted, the
default pipeline is `min` → `eq` → `awh`. Use `stage_overrides`
(see [Imports](#imports)) to change a single stage of an imported or the
default pipeline.

### Global mdrun flags

Machine-specific flags can be set once for all stages under a top-level
`mdrun` key:

```yaml
mdrun:
  bonded: gpu
  ntomp: 8
  ntmpi: 1
```

They are passed to every stage except energy minimizations (stages whose
MDP uses `integrator = steep`, `cg` or `l-bfgs`). A stage's own `mdrun`
flags take precedence. Because these flags describe the machine rather than
the simulation, changing them does not rerun finished steps.

`bonded: gpu` is left out for simulations whose target molecule is a single
bead, because mdrun refuses `-bonded gpu` when there are no bonded
interactions to compute.

## Large molecule libraries

A single ITP holding tens of thousands of molecules is too slow to read for
every simulation. Split it once into a library of per-molecule files:

```bash
martini_harness split-itp all_molecules.itp -o library/
```

This writes `library/index.tsv` (molecule name → file) and one ITP per
molecule, sharded into 256 subdirectories. Splitting 77,000 molecules
(170 MB) takes about 10 s, most of it creating the files. Duplicate names
keep the first definition, content before the first `[ moleculetype ]` is
saved to `_preamble.itp` (not included automatically), and blocks with
unbalanced `#if`/`#endif` are reported.

Pick a random subset into a file that a campaign can import:

```bash
martini_harness sample library/ -n 500 --seed 1 -o campaigns/screen/targets.yaml
```

The file sets `molecule_library` and the `target` axis of the matrix. Without
`--seed` a random seed is chosen and recorded in the file.

For names formatted `<length>_<molid>`, `--per-group K` replaces `-n` and
draws exactly K molecules of every length (the group is the part of the name
before the first `_`). Group sizes are logged, and sampling fails if a group
has fewer than K molecules:

```bash
martini_harness sample library/ --per-group 50 --seed 1 -o campaigns/screen/targets.yaml
``` A campaign then
only adds what is missing:

```yaml
# campaigns/screen/campaign.yaml
imports:
  - targets.yaml
  - ../../examples/shared/solvation_pipeline.yaml
output_dir: out
matrix:
  solvent: [W, OCO]
```

Targets and solvents are looked up in `molecules` first, then in the
topology includes, then in `molecule_library`. Only the index and the ITPs
of the sampled molecules are read when the campaign is built or run, and each
`system.top` includes just its own molecule's file.

## Output and resuming

Each simulation gets a directory `<output_dir>/<name>/` containing
`system.top`, `init.gro`, the rendered `<stage>.mdp`, all GROMACS outputs, the
program output (`<stage>.grompp.out`, `<stage>.mdrun.out`, `coords.log`), the
resolved parameters (`simulation.json`) and the state file
(`.harness_state.json`).

Each step gets a key: a hash of its inputs (topology and included ITP files,
rendered MDP, flags, density, ...) chained with the previous step's key. A
step is skipped when the state file records it as done with the same key and
its outputs still exist; once a step runs, every later step runs too. So
editing a stage's MDP reruns that stage and everything after it, while
changing `density` reruns everything from `coords` on.

A failing simulation is skipped and the campaign continues with the next
one; only Ctrl-C stops the loop. A failed step is recorded in the state file
(`status` shows the error) and is retried on the next `run`. Problems outside
a step, such as an unreadable state file or a simulation directory that
cannot be created, also skip only that simulation.
`run` exits with code 1 if any simulation failed and 2 for configuration
errors.

## Checking production runs

`check` counts how many simulations finished their production stage (by
default the last stage, e.g. `awh`; choose another with `--stage`):

```bash
martini_harness check campaigns/screen/campaign.yaml
```

```
Production stage 'awh': 1000 simulation(s)
  complete        612  ( 61.2%)
  incomplete        4  (  0.4%)  running or interrupted, median progress 37%
  failed           21  (  2.1%)  at awh 15, coords 6
      reasons: LINCS, then non-finite energy 9, GPU 6, coordinate generation (polyply) 6
  not started     363  ( 36.3%)  next step: eq 1, topology 362
```

A run is **complete** when the last frame of its dH/dλ output
(`<stage>.xvg`) has been written: its last time reaches
`tinit + ⌊nsteps / nstdhdl⌋ · nstdhdl · dt`, taken from the
`<stage>_mdout.mdp` grompp wrote (an `-nsteps` mdrun flag is honoured). Only
the end of each xvg file is read, and a half-written last line is ignored.
Runs whose xvg stops early are **incomplete** (still running or
interrupted); **failed** and **not started** come from the state files.

### Failure reasons

When a step fails, the harness reads the end of its program output
(`<stage>.mdrun.out`, `<stage>.log`, `<stage>.grompp.out`), extracts the
GROMACS error message and classifies it, e.g. `LINCS`, `LINCS, then
non-finite energy`, `GPU`, `non-finite forces (pforce)`, `too many grompp
warnings`, `atoms moved too far`, `out of memory`, `killed by signal 9` or
`coordinate generation (polyply)`. LINCS warnings and large forces listed by
`-pforce` before the error are added as notes. The reason is stored in the
state file and shown at the end of `run`, by `status` and by `check`. `check`
also reads the logs directly, so it gives reasons for older batches and for
incomplete runs that crashed without the failure being recorded.

`--details` adds one line per simulation with its progress and reason; `--list
complete|incomplete|failed|not-started` prints just the names, e.g. to
feed into `run --only`.

## Analyses

`analyze` runs a program (a GROMACS tool or anything else) in the directory
of every simulation whose stage finished and collects a result from one of
its output files into a table: the last data line, or the signed difference
between the first and the last line. Analyses are defined in the campaign
file, so they can be imported like the rest of it:

```yaml
analyses:
  - name: dG
    stage: awh                          # default: the last stage
    command: gmx awh -s awh.tpr -f awh.edr -o {tmpdir}/awh
    output: "{tmpdir}/awh_t{end_time}.xvg"
    result: first-last                  # PMF(first λ state) - PMF(last λ state)
    columns: [dG_solv]
```

```bash
martini_harness analyze campaigns/screen/campaign.yaml   # all analyses
martini_harness analyze campaign.yaml --analysis dG --only row42_in_W
```

| key       | meaning                                                                      |
|-----------|------------------------------------------------------------------------------|
| `name`    | Analysis name; results go to `<output_dir>/<name>.tsv`                       |
| `stage`   | Stage to analyse; simulations without `<stage>.gro` are skipped (default: last stage) |
| `command` | Program and arguments, as a string (split like a shell command) or a list    |
| `output`  | File, relative to the simulation directory, that holds the result          |
| `result`  | `last` (default): all columns of the last data line. `first-last` / `last-first`: signed difference of one column between the first and the last data line |
| `column`  | Column the difference is taken of (0 = x; default 1, the first y column)     |
| `stdin`   | Text sent to the program, e.g. a group selection for `gmx energy`            |
| `columns` | Names of the result columns; a difference has one (default: the analysis name) |

`command`, `output` and `stdin` may use the matrix variables (`{target}`,
`{solvent}`, ...) and:

- `{stage}` and `{workdir}` (the simulation directory, where the command runs),
- `{tmpdir}`, a scratch directory inside the simulation directory that is
  removed afterwards, for tools like `gmx awh` that write one file per frame,
- `{end_time}`, the simulated time of the stage in ps, `tinit + nsteps · dt`
  from `<stage>_mdout.mdp`, formatted the way GROMACS names its files
  (`10000`, `1e+06`). With the example above it selects the PMF written at
  the end of the run.

**Free energy from AWH.** `gmx awh` normalises the PMF so that its minimum is
0, so the last value alone is not the free energy. The solvation free energy
is the PMF of the coupled state minus that of the decoupled state. In the
default `awh.mdp`, `vdw_lambdas` runs from 1 (state 0, coupled) to 0 (last
state, decoupled), so `result: first-last` gives ΔG_solv in kJ/mol. For λ
values that run from 0 (decoupled) to 1 (coupled), as in the GROMACS
solvation tutorial, use `last-first`. The difference is formatted with 6
significant digits.

The table has one row per simulation: its name and variables, `status`
(`ok`, `not-finished` or `failed`), the result values and an
`error` column. The program's output is saved to `<name>.analysis.out` in
the simulation directory. `analyze` exits with code 1 if the analysis failed
for any simulation.

## Development

```bash
pytest                    # all tests
pytest -m "not gromacs"   # skip the end-to-end test that needs gmx
```

Most tests use a fake backend (see `tests/conftest.py`) instead of polyply and
GROMACS. `martini_harness/backend.py` is the only module that calls them.

## Package layout

```
martini_harness/
├── cli.py          # command-line interface
├── config.py       # YAML schema (pydantic) and loading
├── matrix.py       # build-matrix expansion
├── simulation.py   # simulations, steps, resume logic, campaign runner
├── check.py        # production completion check (xvg last frame)
├── analysis.py     # analysis programs run per simulation, results table
├── diagnose.py     # failure reasons from GROMACS output
├── state.py        # per-simulation state file
├── backend.py      # polyply / GROMACS execution
├── mdp.py          # MDP rendering and templating
├── topology.py     # topology writing, ITP parsing
├── library.py      # splitting large ITPs, molecule libraries, sampling
└── data/           # Martini 3 ITPs and MDP templates
examples/           # sample campaigns, shared import files, molecule
tests/
```
