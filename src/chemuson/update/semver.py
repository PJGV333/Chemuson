"""Utilidades de versionado semántico para updates."""

from __future__ import annotations

import re
from dataclasses import dataclass
from itertools import zip_longest
from typing import Iterable

from chemuson.update.types import UpdateChannel, coerce_update_channel

_SEMVER_RE = re.compile(
    r"^(?P<prefix>v)?"
    r"(?P<major>0|[1-9]\d*)\."
    r"(?P<minor>0|[1-9]\d*)\."
    r"(?P<patch>0|[1-9]\d*)"
    r"(?:-(?P<prerelease>[0-9A-Za-z.-]+))?"
    r"(?:\+(?P<build>[0-9A-Za-z.-]+))?$"
)


@dataclass(frozen=True, slots=True)
class SemVer:
    """Representación canónica de una versión SemVer."""

    major: int
    minor: int
    patch: int
    prerelease: tuple[str, ...] = ()
    build: str = ""

    def normalized(self) -> str:
        """Devuelve representación normalizada sin prefijo `v`."""
        base = f"{self.major}.{self.minor}.{self.patch}"
        if self.prerelease:
            base += "-" + ".".join(self.prerelease)
        if self.build:
            base += "+" + self.build
        return base

    @property
    def is_prerelease(self) -> bool:
        """Indica si corresponde a prerelease."""
        return bool(self.prerelease)


def parse_semver(value: str) -> SemVer:
    """Parsea una versión semántica (admite prefijo `v`)."""
    raw = str(value or "").strip()
    match = _SEMVER_RE.match(raw)
    if not match:
        raise ValueError(f"Versión SemVer inválida: {value!r}")
    prerelease_raw = match.group("prerelease") or ""
    prerelease = tuple(part for part in prerelease_raw.split(".") if part) if prerelease_raw else ()
    return SemVer(
        major=int(match.group("major")),
        minor=int(match.group("minor")),
        patch=int(match.group("patch")),
        prerelease=prerelease,
        build=match.group("build") or "",
    )


def _is_numeric_identifier(part: str) -> bool:
    return part.isdigit()


def _compare_prerelease(a: Iterable[str], b: Iterable[str]) -> int:
    """Compara prereleases según SemVer 2.0.0."""
    a_parts = list(a)
    b_parts = list(b)
    if not a_parts and not b_parts:
        return 0
    if not a_parts:
        return 1
    if not b_parts:
        return -1
    for left, right in zip_longest(a_parts, b_parts, fillvalue=None):
        if left is None:
            return -1
        if right is None:
            return 1
        left_num = _is_numeric_identifier(left)
        right_num = _is_numeric_identifier(right)
        if left_num and right_num:
            li = int(left)
            ri = int(right)
            if li < ri:
                return -1
            if li > ri:
                return 1
            continue
        if left_num and not right_num:
            return -1
        if right_num and not left_num:
            return 1
        if left < right:
            return -1
        if left > right:
            return 1
    return 0


def compare_versions(left: str, right: str) -> int:
    """Compara dos versiones semánticas y devuelve -1, 0, 1."""
    a = parse_semver(left)
    b = parse_semver(right)
    for av, bv in ((a.major, b.major), (a.minor, b.minor), (a.patch, b.patch)):
        if av < bv:
            return -1
        if av > bv:
            return 1
    return _compare_prerelease(a.prerelease, b.prerelease)


def is_newer_version(candidate: str, current: str) -> bool:
    """Indica si `candidate` es más reciente que `current`."""
    return compare_versions(candidate, current) > 0


def is_prerelease(version: str) -> bool:
    """Indica si una versión es prerelease."""
    return parse_semver(version).is_prerelease


def channel_accepts_version(channel: UpdateChannel | str, version: str) -> bool:
    """Valida si el canal permite consumir la versión indicada."""
    ch = coerce_update_channel(channel)
    if ch == UpdateChannel.STABLE:
        return not is_prerelease(version)
    return True
