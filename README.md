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
martini_harness list-molecules examples/row42.itp             # molecule types in ITPs
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

A failing step stops its simulation and is recorded in the state file
(`status` shows the error); the campaign continues with the next simulation.
`run` exits with code 1 if any simulation failed and 2 for configuration
errors.

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
├── state.py        # per-simulation state file
├── backend.py      # polyply / GROMACS execution
├── mdp.py          # MDP rendering and templating
├── topology.py     # topology writing, ITP parsing
└── data/           # Martini 3 ITPs and MDP templates
examples/           # sample campaigns, shared import files, molecule
tests/
```
