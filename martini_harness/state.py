"""Per-simulation state file used to resume interrupted campaigns."""

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

STATE_FILE_NAME = ".harness_state.json"
STATE_VERSION = 1

DONE = "done"
FAILED = "failed"


class StateFile:
    """Records which steps of a simulation finished, keyed by their inputs.

    Each step entry stores the step's cache key (a hash of everything that
    determines its output, chained with the previous step's key), so a step
    counts as done only if it finished with exactly the current inputs.
    """

    def __init__(self, workdir: Path):
        self.path = Path(workdir) / STATE_FILE_NAME
        self._steps: dict[str, dict[str, Any]] = {}
        if self.path.exists():
            data = json.loads(self.path.read_text())
            if data.get("version") == STATE_VERSION:
                self._steps = data.get("steps", {})

    def get(self, step: str) -> dict[str, Any] | None:
        return self._steps.get(step)

    def is_done(self, step: str, key: str) -> bool:
        entry = self._steps.get(step)
        return entry is not None and entry["status"] == DONE and entry["key"] == key

    def record(self, step: str, key: str, status: str, error: str | None = None) -> None:
        entry = {
            "key": key,
            "status": status,
            "time": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        if error is not None:
            entry["error"] = error
        self._steps[step] = entry
        self._save()

    def discard(self, steps: list[str]) -> None:
        for step in steps:
            self._steps.pop(step, None)
        self._save()

    def _save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"version": STATE_VERSION, "steps": self._steps}, indent=2))
        os.replace(tmp, self.path)
