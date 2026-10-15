"""Execution backends that run the external programs (polyply, GROMACS)."""

import logging
from pathlib import Path
from typing import Any, Protocol


class BackendError(RuntimeError):
    """Raised when an external program fails."""


class Backend(Protocol):
    """Runs the external programs of a simulation.

    All paths passed in are absolute; ``log`` receives the program's
    stdout/stderr.
    """

    def generate_coordinates(
        self, topology: Path, output: Path, density: float, name: str, log: Path
    ) -> None: ...

    def grompp(self, flags: dict[str, Any], log: Path) -> None: ...

    def mdrun(self, flags: dict[str, Any], log: Path) -> None: ...


class GromacsWrapperBackend:
    """Runs polyply in-process and GROMACS through GromacsWrapper."""

    def generate_coordinates(
        self, topology: Path, output: Path, density: float, name: str, log: Path
    ) -> None:
        from polyply import gen_coords

        handler = logging.FileHandler(log, mode="w")
        root = logging.getLogger()
        root.addHandler(handler)
        try:
            gen_coords(toppath=topology, outpath=output, name=name, density=density)
        except Exception as exc:
            raise BackendError(f"polyply gen_coords failed: {exc}; see {log}") from exc
        finally:
            root.removeHandler(handler)
            handler.close()

    def grompp(self, flags: dict[str, Any], log: Path) -> None:
        import gromacs

        self._run(gromacs.grompp, flags, log)

    def mdrun(self, flags: dict[str, Any], log: Path) -> None:
        import gromacs

        self._run(gromacs.mdrun, flags, log)

    @staticmethod
    def _run(tool: Any, flags: dict[str, Any], log: Path) -> None:
        flags = {key: str(value) if isinstance(value, Path) else value for key, value in flags.items()}
        with open(log, "w") as handle:
            try:
                returncode, _, _ = tool(stdout=handle, stderr=handle, **flags)
            except Exception as exc:
                raise BackendError(f"{tool.command_name} failed; see {log}") from exc
        if returncode != 0:
            raise BackendError(f"{tool.command_name} exited with code {returncode}; see {log}")
