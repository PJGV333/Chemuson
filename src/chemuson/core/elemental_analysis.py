"""Elemental analysis calculations for formula strings.

This module intentionally has no GUI dependencies.  It parses molecular
formula/adduct strings, calculates average molecular weight and elemental
percentages, compares calculated vs found percentages, and ranks transparent
hydrate/solvate candidates.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from itertools import combinations, product
import math
from typing import Iterable, Mapping, Sequence


DEFAULT_TOLERANCE = 0.40
MIDDLE_DOT = "\u00b7"
_EPSILON = 1e-10


class FormulaError(ValueError):
    """Raised when a formula cannot be parsed into element counts."""


# Average atomic weights.  Values are intentionally local and deterministic.
# They are sufficient for common organic/inorganic elemental analysis work.
ATOMIC_WEIGHTS: dict[str, float] = {
    "H": 1.00794,
    "He": 4.002602,
    "Li": 6.941,
    "Be": 9.012182,
    "B": 10.811,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "F": 18.9984032,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.3050,
    "Al": 26.9815386,
    "Si": 28.0855,
    "P": 30.973762,
    "S": 32.065,
    "Cl": 35.453,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955912,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938045,
    "Fe": 55.845,
    "Co": 58.933195,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.64,
    "As": 74.92160,
    "Se": 78.96,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90585,
    "Zr": 91.224,
    "Nb": 92.90638,
    "Mo": 95.96,
    "Tc": 98.0,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.411,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.9054519,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pr": 140.90765,
    "Nd": 144.242,
    "Pm": 145.0,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.92535,
    "Dy": 162.500,
    "Ho": 164.93032,
    "Er": 167.259,
    "Tm": 168.93421,
    "Yb": 173.054,
    "Lu": 174.9668,
    "Hf": 178.49,
    "Ta": 180.94788,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.966569,
    "Hg": 200.59,
    "Tl": 204.3833,
    "Pb": 207.2,
    "Bi": 208.98040,
    "Po": 209.0,
    "At": 210.0,
    "Rn": 222.0,
    "Fr": 223.0,
    "Ra": 226.0,
    "Ac": 227.0,
    "Th": 232.03806,
    "Pa": 231.03588,
    "U": 238.02891,
    "Np": 237.0,
    "Pu": 244.0,
    "Am": 243.0,
    "Cm": 247.0,
    "Bk": 247.0,
    "Cf": 251.0,
    "Es": 252.0,
    "Fm": 257.0,
    "Md": 258.0,
    "No": 259.0,
    "Lr": 262.0,
    "Rf": 267.0,
    "Db": 268.0,
    "Sg": 271.0,
    "Bh": 272.0,
    "Hs": 270.0,
    "Mt": 276.0,
    "Ds": 281.0,
    "Rg": 280.0,
    "Cn": 285.0,
    "Nh": 284.0,
    "Fl": 289.0,
    "Mc": 288.0,
    "Lv": 293.0,
    "Ts": 294.0,
    "Og": 294.0,
}

ELEMENT_SYMBOLS = frozenset(ATOMIC_WEIGHTS)


@dataclass(frozen=True)
class SolventEntry:
    """Named solvent/adduct entry."""

    name: str
    formula: str


SOLVENT_LIBRARY: dict[str, SolventEntry] = {
    "H2O": SolventEntry("H2O", "H2O"),
    "MeOH": SolventEntry("MeOH", "CH4O"),
    "EtOH": SolventEntry("EtOH", "C2H6O"),
    "acetone": SolventEntry("acetone", "C3H6O"),
    "MeCN": SolventEntry("MeCN", "C2H3N"),
    "DMSO": SolventEntry("DMSO", "C2H6OS"),
    "DMF": SolventEntry("DMF", "C3H7NO"),
    "THF": SolventEntry("THF", "C4H8O"),
    "CH2Cl2": SolventEntry("CH2Cl2", "CH2Cl2"),
    "CHCl3": SolventEntry("CHCl3", "CHCl3"),
    "EtOAc": SolventEntry("EtOAc", "C4H8O2"),
    "hexane": SolventEntry("hexane", "C6H14"),
    "toluene": SolventEntry("toluene", "C7H8"),
    "diethyl ether": SolventEntry("diethyl ether", "C4H10O"),
}


@dataclass(frozen=True)
class ComparisonEntry:
    """One calculated/found comparison row."""

    element: str
    calculated: float
    found: float | None
    delta: float | None
    tolerance: float
    passed: bool | None


@dataclass(frozen=True)
class SolvateCandidate:
    """Ranked candidate hydrate/solvate formula."""

    formula: str
    molecular_weight: float
    calculated: dict[str, float]
    deltas: dict[str, float]
    score: float
    max_absolute_delta: float
    passed: bool
    solvent_equivalents: dict[str, float]
    warning: str = ""


def parse_formula(formula: str) -> dict[str, float]:
    """Parse a molecular formula/adduct expression into element counts.

    Supports parentheses, integer/decimal stoichiometric coefficients, middle
    dot or period adduct notation, and plus-separated adducts with optional
    leading equivalents, for example ``CuSO4\u00b75H2O`` or
    ``C20H18N4O2 + 0.25 CH3OH``.
    """

    if formula is None or not str(formula).strip():
        raise FormulaError("Formula is empty.")
    expression = _normalize_formula_text(str(formula))
    components = _split_adduct_components(expression)
    if not components:
        raise FormulaError("Formula is empty.")

    total: dict[str, float] = {}
    for component in components:
        factor, body = _leading_factor(component)
        parser = _FormulaParser(body)
        parsed = parser.parse()
        _merge_counts(total, parsed, factor)
    return _clean_composition(total)


def molecular_weight(composition: Mapping[str, float]) -> float:
    """Return average molecular weight for an element composition."""

    total = 0.0
    for element, raw_count in composition.items():
        count = _validated_count(element, raw_count)
        if count <= _EPSILON:
            continue
        try:
            weight = ATOMIC_WEIGHTS[element]
        except KeyError as exc:
            raise FormulaError(f"Atomic weight is not available for '{element}'.") from exc
        total += weight * count
    if total <= _EPSILON:
        raise FormulaError("Composition has no positive atom counts.")
    return total


def elemental_percentages(composition: Mapping[str, float]) -> dict[str, float]:
    """Return elemental mass percentages for an element composition."""

    mw = molecular_weight(composition)
    percentages: dict[str, float] = {}
    for element, count in composition.items():
        if count <= _EPSILON:
            continue
        percentages[element] = (ATOMIC_WEIGHTS[element] * float(count) / mw) * 100.0
    return percentages


def compare_found_calculated(
    calculated: Mapping[str, float],
    found: Mapping[str, float | None],
    tolerances: Mapping[str, float] | float | None = None,
) -> list[ComparisonEntry]:
    """Compare calculated and found percentages element-by-element."""

    if found:
        elements = list(found.keys())
    else:
        elements = list(calculated.keys())
    rows: list[ComparisonEntry] = []
    for element in _element_order(elements):
        calc_value = float(calculated.get(element, 0.0))
        found_value = _optional_float(found.get(element)) if found else None
        tolerance = _tolerance_for(element, tolerances)
        if found_value is None:
            rows.append(
                ComparisonEntry(
                    element=element,
                    calculated=calc_value,
                    found=None,
                    delta=None,
                    tolerance=tolerance,
                    passed=None,
                )
            )
            continue
        delta = found_value - calc_value
        rows.append(
            ComparisonEntry(
                element=element,
                calculated=calc_value,
                found=found_value,
                delta=delta,
                tolerance=tolerance,
                passed=abs(delta) <= tolerance + _EPSILON,
            )
        )
    return rows


def format_publication_line(
    formula: str,
    calculated: Mapping[str, float],
    found: Mapping[str, float | None],
    selected_elements: Sequence[str],
    decimals: int = 2,
) -> str:
    """Return a manuscript-ready found/calculated elemental analysis line."""

    decimals = max(0, int(decimals))
    elements = [element for element in selected_elements if element in calculated or element in found]
    if not elements:
        elements = _element_order(calculated.keys())
    calc_parts = [
        f"{element}, {_format_number(float(calculated.get(element, 0.0)), decimals)}"
        for element in elements
    ]
    found_parts = []
    for element in elements:
        value = _optional_float(found.get(element))
        rendered = "N/D" if value is None else _format_number(value, decimals)
        found_parts.append(f"{element}, {rendered}")
    return (
        f"Anal. Calcd for {formula}: "
        f"{'; '.join(calc_parts)}. Found: {'; '.join(found_parts)}."
    )


def format_formula(composition: Mapping[str, float], decimals: int = 4) -> str:
    """Format a normalized formula using Hill ordering."""

    parts: list[str] = []
    for element in _element_order(composition.keys()):
        count = float(composition[element])
        if count <= _EPSILON:
            continue
        parts.append(element + _format_count(count, decimals=decimals))
    return "".join(parts)


def find_solvate_candidates(
    base_formula: str,
    found_percentages: Mapping[str, float | None],
    candidate_solvents: Sequence[str | SolventEntry | tuple[str, str]],
    equivalent_range: tuple[float, float] = (0.25, 2.0),
    step_size: float = 0.25,
    selected_elements: Sequence[str] | None = None,
    tolerances: Mapping[str, float] | float | None = None,
    max_components: int = 1,
    max_results: int | None = None,
) -> list[SolvateCandidate]:
    """Grid-search ranked hydrate/solvate candidate formulas.

    Returned formulas are candidate formulas, not definitive assignments.
    """

    base_composition = parse_formula(base_formula)
    base_label = format_formula(base_composition)
    elements = list(selected_elements or found_percentages.keys())
    elements = [element for element in _element_order(elements) if _optional_float(found_percentages.get(element)) is not None]
    if not elements:
        raise FormulaError("At least one selected element needs a found percentage.")

    solvent_entries = [_resolve_solvent_entry(item) for item in candidate_solvents]
    solvent_entries = _dedupe_solvents(solvent_entries)
    max_components = max(0, min(int(max_components), len(solvent_entries)))
    equivalents = _equivalent_values(equivalent_range, step_size)

    candidates: list[SolvateCandidate] = []
    seen_formulas: set[str] = set()

    def add_candidate(parts: Sequence[tuple[SolventEntry, float]]) -> None:
        composition = dict(base_composition)
        formula_parts = [base_label]
        solvent_equivalents: dict[str, float] = {}
        for solvent, equivalent in parts:
            if equivalent <= _EPSILON:
                continue
            solvent_composition = parse_formula(solvent.formula)
            _merge_counts(composition, solvent_composition, equivalent)
            formula_parts.append(f"{_format_count(equivalent, decimals=4)}{solvent.formula}")
            solvent_equivalents[solvent.name] = equivalent
        formula = MIDDLE_DOT.join(formula_parts)
        if formula in seen_formulas:
            return
        seen_formulas.add(formula)
        calculated = elemental_percentages(composition)
        deltas = {
            element: float(found_percentages[element]) - calculated.get(element, 0.0)
            for element in elements
        }
        tolerance_values = {element: _tolerance_for(element, tolerances) for element in elements}
        score = math.sqrt(
            sum((deltas[element] / tolerance_values[element]) ** 2 for element in elements)
            / len(elements)
        )
        max_absolute_delta = max(abs(value) for value in deltas.values())
        passed = all(abs(deltas[element]) <= tolerance_values[element] + _EPSILON for element in elements)
        warning = _candidate_warning(len(solvent_equivalents), len(elements))
        candidates.append(
            SolvateCandidate(
                formula=formula,
                molecular_weight=molecular_weight(composition),
                calculated=calculated,
                deltas=deltas,
                score=score,
                max_absolute_delta=max_absolute_delta,
                passed=passed,
                solvent_equivalents=solvent_equivalents,
                warning=warning,
            )
        )

    add_candidate(())
    for component_count in range(1, max_components + 1):
        for solvent_combo in combinations(solvent_entries, component_count):
            for equivalent_combo in product(equivalents, repeat=component_count):
                if all(value <= _EPSILON for value in equivalent_combo):
                    continue
                add_candidate(tuple(zip(solvent_combo, equivalent_combo)))

    candidates.sort(key=lambda item: (item.score, item.max_absolute_delta, item.molecular_weight, item.formula))
    if max_results is not None:
        return candidates[: max(0, int(max_results))]
    return candidates


def _normalize_formula_text(formula: str) -> str:
    normalized = formula.strip()
    normalized = normalized.replace("[", "(").replace("]", ")")
    normalized = normalized.replace("{", "(").replace("}", ")")
    return normalized


def _split_adduct_components(expression: str) -> list[str]:
    components: list[str] = []
    depth = 0
    start = 0
    i = 0
    while i < len(expression):
        char = expression[i]
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise FormulaError("Unmatched closing parenthesis.")
        elif depth == 0 and (char == "+" or char == MIDDLE_DOT or _is_period_separator(expression, i)):
            component = expression[start:i].strip()
            if not component:
                raise FormulaError("Missing formula component in adduct expression.")
            components.append(component)
            start = i + 1
        i += 1
    if depth != 0:
        raise FormulaError("Unmatched opening parenthesis.")
    component = expression[start:].strip()
    if not component:
        raise FormulaError("Missing formula component in adduct expression.")
    components.append(component)
    return components


def _is_period_separator(expression: str, index: int) -> bool:
    if expression[index] != ".":
        return False
    previous_char = expression[index - 1] if index > 0 else ""
    next_char = expression[index + 1] if index + 1 < len(expression) else ""
    return not (previous_char.isdigit() and next_char.isdigit())


def _leading_factor(component: str) -> tuple[float, str]:
    index = 0
    while index < len(component) and component[index].isspace():
        index += 1
    end = _scan_number(component, index)
    if end == index:
        return 1.0, component.strip()
    factor = _parse_positive_number(component[index:end], context="leading coefficient")
    body = component[end:].strip()
    if not body:
        raise FormulaError("Leading coefficient is missing a formula.")
    return factor, body


class _FormulaParser:
    def __init__(self, text: str) -> None:
        self.text = text.replace(" ", "")
        self.index = 0

    def parse(self) -> dict[str, float]:
        if not self.text:
            raise FormulaError("Formula component is empty.")
        counts = self._parse_group(stop_char=None)
        if self.index != len(self.text):
            raise FormulaError(f"Unexpected text at position {self.index + 1}.")
        return counts

    def _parse_group(self, stop_char: str | None) -> dict[str, float]:
        counts: dict[str, float] = {}
        while self.index < len(self.text):
            char = self.text[self.index]
            if stop_char and char == stop_char:
                return counts
            if char == "(":
                self.index += 1
                inner = self._parse_group(stop_char=")")
                if self.index >= len(self.text) or self.text[self.index] != ")":
                    raise FormulaError("Unmatched opening parenthesis.")
                self.index += 1
                multiplier = self._parse_optional_number(default=1.0)
                _merge_counts(counts, inner, multiplier)
                continue
            if char == ")":
                if stop_char is None:
                    raise FormulaError("Unmatched closing parenthesis.")
                return counts
            if char.isupper():
                element = self._parse_element()
                multiplier = self._parse_optional_number(default=1.0)
                _merge_counts(counts, {element: 1.0}, multiplier)
                continue
            if char.isdigit() or char == ".":
                raise FormulaError(
                    f"Unexpected coefficient at position {self.index + 1}; "
                    "coefficients must follow an element or group."
                )
            raise FormulaError(f"Unexpected character '{char}' at position {self.index + 1}.")
        if stop_char:
            raise FormulaError("Unmatched opening parenthesis.")
        return counts

    def _parse_element(self) -> str:
        start = self.index
        self.index += 1
        if self.index < len(self.text) and self.text[self.index].islower():
            self.index += 1
        element = self.text[start:self.index]
        if element not in ELEMENT_SYMBOLS:
            raise FormulaError(f"Unknown element symbol '{element}'.")
        return element

    def _parse_optional_number(self, default: float) -> float:
        start = self.index
        end = _scan_number(self.text, start)
        if end == start:
            return default
        self.index = end
        return _parse_positive_number(self.text[start:end], context="stoichiometric coefficient")


def _scan_number(text: str, start: int) -> int:
    index = start
    saw_digit = False
    saw_decimal = False
    while index < len(text):
        char = text[index]
        if char.isdigit():
            saw_digit = True
            index += 1
            continue
        if char == "." and not saw_decimal:
            saw_decimal = True
            index += 1
            continue
        break
    if not saw_digit:
        return start
    return index


def _parse_positive_number(text: str, context: str) -> float:
    try:
        value = float(Decimal(text))
    except (InvalidOperation, ValueError) as exc:
        raise FormulaError(f"Invalid {context} '{text}'.") from exc
    if not math.isfinite(value) or value <= 0:
        raise FormulaError(f"Invalid {context} '{text}'; value must be positive.")
    return value


def _merge_counts(target: dict[str, float], source: Mapping[str, float], factor: float = 1.0) -> None:
    for element, count in source.items():
        target[element] = target.get(element, 0.0) + float(count) * float(factor)


def _clean_composition(composition: Mapping[str, float]) -> dict[str, float]:
    cleaned: dict[str, float] = {}
    for element, count in composition.items():
        value = _validated_count(element, count)
        if value <= _EPSILON:
            continue
        rounded = round(value)
        cleaned[element] = float(rounded) if abs(value - rounded) <= _EPSILON else value
    if not cleaned:
        raise FormulaError("Formula has no positive atom counts.")
    return dict(sorted(cleaned.items(), key=lambda item: _element_sort_key(item[0])))


def _validated_count(element: str, raw_count: float) -> float:
    if element not in ELEMENT_SYMBOLS:
        raise FormulaError(f"Unknown element symbol '{element}'.")
    try:
        count = float(raw_count)
    except (TypeError, ValueError) as exc:
        raise FormulaError(f"Invalid atom count for '{element}'.") from exc
    if not math.isfinite(count) or count < 0:
        raise FormulaError(f"Invalid atom count for '{element}'; value must be non-negative.")
    return count


def _element_order(elements: Iterable[str]) -> list[str]:
    return sorted(dict.fromkeys(elements), key=_element_sort_key)


def _element_sort_key(element: str) -> tuple[int, str]:
    if element == "C":
        return (0, element)
    if element == "H":
        return (1, element)
    return (2, element)


def _format_count(value: float, decimals: int) -> str:
    if abs(value - 1.0) <= _EPSILON:
        return ""
    rounded = round(value)
    if abs(value - rounded) <= _EPSILON:
        return str(int(rounded))
    return _format_number(value, decimals).rstrip("0").rstrip(".")


def _format_number(value: float, decimals: int) -> str:
    return f"{value:.{max(0, int(decimals))}f}"


def _optional_float(value: float | None) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _tolerance_for(element: str, tolerances: Mapping[str, float] | float | None) -> float:
    if tolerances is None:
        value = DEFAULT_TOLERANCE
    elif isinstance(tolerances, Mapping):
        value = float(tolerances.get(element, DEFAULT_TOLERANCE))
    else:
        value = float(tolerances)
    if not math.isfinite(value) or value <= 0:
        raise FormulaError(f"Tolerance for '{element}' must be positive.")
    return value


def _resolve_solvent_entry(item: str | SolventEntry | tuple[str, str]) -> SolventEntry:
    if isinstance(item, SolventEntry):
        return item
    if isinstance(item, tuple):
        name, formula = item
        parse_formula(formula)
        return SolventEntry(str(name).strip() or str(formula), str(formula).strip())
    key = str(item).strip()
    if key in SOLVENT_LIBRARY:
        return SOLVENT_LIBRARY[key]
    parse_formula(key)
    return SolventEntry(key, key)


def _dedupe_solvents(solvents: Sequence[SolventEntry]) -> list[SolventEntry]:
    seen: set[tuple[str, str]] = set()
    result: list[SolventEntry] = []
    for solvent in solvents:
        key = (solvent.name, solvent.formula)
        if key in seen:
            continue
        seen.add(key)
        result.append(solvent)
    return result


def _equivalent_values(equivalent_range: tuple[float, float], step_size: float) -> list[float]:
    if len(equivalent_range) != 2:
        raise FormulaError("Equivalent range must contain minimum and maximum values.")
    start = Decimal(str(equivalent_range[0]))
    stop = Decimal(str(equivalent_range[1]))
    step = Decimal(str(step_size))
    if start < 0 or stop < 0 or step <= 0:
        raise FormulaError("Equivalent range and step must be non-negative with a positive step.")
    if stop < start:
        raise FormulaError("Equivalent range maximum must be greater than or equal to minimum.")
    values: list[float] = []
    current = start
    guard = 0
    while current <= stop + Decimal("1e-12"):
        values.append(float(current))
        current += step
        guard += 1
        if guard > 10000:
            raise FormulaError("Equivalent grid is too large.")
    return values


def _candidate_warning(component_count: int, element_count: int) -> str:
    warnings: list[str] = []
    if component_count > 1:
        warnings.append("Multiple solvent components can overfit elemental percentages.")
    if component_count >= max(1, element_count - 1):
        warnings.append("Fitted variables approach the number of experimental elements.")
    return " ".join(warnings)
