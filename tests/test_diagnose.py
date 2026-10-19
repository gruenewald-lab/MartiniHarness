import os
import time
from pathlib import Path

import pytest

from martini_harness.diagnose import diagnose_step, diagnose_text, gromacs_error

# Tails of real GROMACS 2024.2 output, see tests/data/gromacs_errors.
SAMPLES = Path(__file__).parent / "data" / "gromacs_errors"


@pytest.mark.parametrize(
    "sample, reason, message_start",
    [
        ("gpu.mdrun.out", "GPU", "Nonbonded interactions on the GPU were requested with -nb gpu"),
        ("bgpu.mdrun.out", "GPU", "Bonded interactions on GPU are not supported"),
        ("lpf.mdrun.out", "non-finite forces (pforce)", "At step 2 detected non-finite forces on 2 atoms"),
        ("lincs.mdrun.out", "LINCS, then non-finite energy", "Step 100: The total potential energy is nan"),
        ("warn.grompp.out", "too many grompp warnings", "Too many warnings (2)."),
    ],
)
def test_real_gromacs_errors(sample, reason, message_start):
    diagnosis = diagnose_text((SAMPLES / sample).read_text())
    assert diagnosis.reason == reason
    assert diagnosis.message.startswith(message_start)


def test_notes_count_lincs_warnings_and_large_forces():
    text = (
        "Step 1, time 0.02 (ps)  LINCS WARNING\n" * 3
        + "step 1 atom     32  x    0.564    0.001    2.289  force  5.52496e+10\n" * 2
        + (SAMPLES / "gpu.mdrun.out").read_text()
    )
    diagnosis = diagnose_text(text)
    assert diagnosis.notes == ["after 3 LINCS warning(s)", "2 large force(s) listed by -pforce"]
    assert str(diagnosis).startswith("GPU (after 3 LINCS warning(s); 2 large force(s)")


def test_explicit_lincs_error_and_interrupted_runs():
    block = (
        "Program:     gmx mdrun\nSource file: constr.cpp\n\nFatal error:\n"
        "Too many LINCS warnings (1000)\nIf you know what you are doing...\n\n"
        "For more information and tips for troubleshooting, please check\n"
    )
    assert diagnose_text(block).reason == "LINCS"
    assert gromacs_error(block).startswith("Too many LINCS warnings (1000) If you know")
    interrupted = diagnose_text("Step 5, time 0.1 (ps)  LINCS WARNING\n")
    assert interrupted.reason == "LINCS" and "interrupted" in interrupted.message
    assert diagnose_text("Segmentation fault (core dumped)\n").reason.startswith("crashed")
    assert diagnose_text("step 0 atom 1 x 1 2 3 force 1e5\nFinished mdrun\n") is None


def test_diagnose_step_uses_current_attempt(tmp_path):
    (tmp_path / "eq.mdrun.out").write_text((SAMPLES / "gpu.mdrun.out").read_text())
    old = time.time() - 100
    os.utime(tmp_path / "eq.mdrun.out", (old, old))  # left over from an earlier attempt
    (tmp_path / "eq.grompp.out").write_text((SAMPLES / "warn.grompp.out").read_text())
    assert diagnose_step(tmp_path, "eq").reason == "too many grompp warnings"

    (tmp_path / "eq.mdrun.out").write_text((SAMPLES / "lpf.mdrun.out").read_text())  # newer
    assert diagnose_step(tmp_path, "eq").reason == "non-finite forces (pforce)"


def test_diagnose_step_signals_and_coordinates(tmp_path):
    assert diagnose_step(tmp_path, "awh", "mdrun exited with code -9; see x").reason.startswith(
        "killed by signal 9"
    )
    coords = diagnose_step(tmp_path, "coords", "polyply gen_coords failed: no space")
    assert coords.reason == "coordinate generation (polyply)" and "no space" in coords.message
    assert diagnose_step(tmp_path, "awh") is None
