"""Rendering of MDP templates and string templating of config values."""

import string
from typing import Any

MINIMIZATION_INTEGRATORS = {"steep", "cg", "l-bfgs"}


class TemplateError(ValueError):
    """Raised when a template references an unknown variable or placeholder."""


def render_value(value: Any, variables: dict[str, Any]) -> Any:
    """Fill ``{variable}`` fields in strings, recursing into lists and dicts.

    Non-string values are returned unchanged.
    """
    if isinstance(value, str):
        fields = [field for _, field, _, _ in string.Formatter().parse(value) if field]
        missing = sorted({field for field in fields if field not in variables})
        if missing:
            raise TemplateError(
                f"unknown variable(s) {', '.join(missing)} in {value!r}; "
                f"available: {', '.join(sorted(variables))}"
            )
        return value.format_map(variables)
    if isinstance(value, dict):
        return {key: render_value(item, variables) for key, item in value.items()}
    if isinstance(value, list):
        return [render_value(item, variables) for item in value]
    return value


def render_mdp(
    template: str,
    substitute: dict[str, str] | None = None,
    params: dict[str, Any] | None = None,
) -> str:
    """Render an MDP file from a template.

    Args:
        template: Content of the template MDP file.
        substitute: Literal placeholder replacements, e.g. ``{"TOCOUPLE": "row42"}``.
            Every placeholder must occur in the template.
        params: MDP parameters to set. Existing entries are replaced (``-`` and
            ``_`` are treated as equivalent, as GROMACS does), others are appended.

    Returns:
        The rendered MDP content.
    """
    text = template
    for placeholder, replacement in (substitute or {}).items():
        if placeholder not in text:
            raise TemplateError(f"placeholder {placeholder!r} not found in MDP template")
        text = text.replace(placeholder, str(replacement))
    return set_mdp_params(text, params or {})


def set_mdp_params(text: str, params: dict[str, Any]) -> str:
    """Set MDP parameters in ``text``, replacing existing entries or appending them."""
    wanted = {_normalize(key): (key, value) for key, value in params.items()}
    found = set()
    lines = []
    for line in text.splitlines():
        body = line.split(";", 1)[0]
        if "=" in body:
            key = _normalize(body.split("=", 1)[0])
            if key in wanted:
                name, value = wanted[key]
                lines.append(_format_param(name, value))
                found.add(key)
                continue
        lines.append(line)
    for key, (name, value) in wanted.items():
        if key not in found:
            lines.append(_format_param(name, value))
    return "\n".join(lines) + "\n"


def get_mdp_param(text: str, key: str) -> str | None:
    """Return the value of ``key`` in MDP ``text`` (last entry wins), or None."""
    wanted = _normalize(key)
    value = None
    for line in text.splitlines():
        name, sep, raw_value = line.split(";", 1)[0].partition("=")
        if sep and _normalize(name) == wanted:
            value = raw_value.strip()
    return value


def is_minimization(text: str) -> bool:
    """Whether the MDP uses an energy-minimization integrator."""
    integrator = get_mdp_param(text, "integrator") or "md"
    return integrator.lower() in MINIMIZATION_INTEGRATORS


def _normalize(key: str) -> str:
    return key.strip().lower().replace("_", "-")


def _format_param(name: str, value: Any) -> str:
    if isinstance(value, bool):
        value = "yes" if value else "no"
    elif isinstance(value, (list, tuple)):
        value = " ".join(str(item) for item in value)
    return f"{name:<24} = {value}"
