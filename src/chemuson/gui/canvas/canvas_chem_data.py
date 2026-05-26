"""Datos químicos estables usados por los mixins del canvas."""

from chemuson.core.model import SIMPLE_HYDROGEN_GROUP_LABELS as CORE_SIMPLE_HYDROGEN_GROUP_LABELS

FUNCTIONAL_GROUP_LABELS = [
    "NH2",
    "NO2",
    "OH",
    "COOH",
    "CO2H",
    "CHO",
    "CN",
    "SO3H",
    "SH",
    "OMe",
    "OEt",
    "Me",
    "Et",
    "iPr",
    "tBu",
    "TBS",
    "Si",
    "Ph",
    "Ac",
    "Boc",
    "Ts",
    "R",
]

FUNCTIONAL_GROUP_ALIASES = {label.lower(): label for label in FUNCTIONAL_GROUP_LABELS}

SIMPLE_HYDROGEN_GROUP_LABELS = dict(CORE_SIMPLE_HYDROGEN_GROUP_LABELS)

CANONICAL_SIMPLE_HYDROGEN_GROUP_LABELS: dict[tuple[str, int], str] = {
    ("N", 2): "NH2",
    ("N", 1): "NH",
    ("O", 1): "OH",
    ("S", 1): "SH",
}

ANALYSIS_MARGIN_PX = 14.0
ANALYSIS_MIN_PEAK_PERCENT = 1.0
ANALYSIS_MAX_PEAKS = 10
ANALYSIS_DIST_KEEP = 400

ATOMIC_WEIGHTS = {
    "H": 1.00794,
    "C": 12.0107,
    "N": 14.0067,
    "O": 15.9994,
    "S": 32.065,
    "P": 30.973762,
    "F": 18.998403,
    "Cl": 35.453,
    "Br": 79.904,
    "I": 126.90447,
    "Si": 28.0855,
    "B": 10.811,
    "Se": 78.96,
    "Li": 6.941,
    "Na": 22.98977,
    "K": 39.0983,
    "Mg": 24.305,
}

MONOISOTOPIC_MASSES = {
    "H": 1.007825032,
    "C": 12.0,
    "N": 14.003074005,
    "O": 15.99491462,
    "S": 31.972071,
    "P": 30.973761998,
    "F": 18.998403163,
    "Cl": 34.96885268,
    "Br": 78.9183376,
    "I": 126.904468,
    "Si": 27.9769265,
    "B": 11.009305,
    "Se": 79.916521,
    "Li": 7.016003,
    "Na": 22.98976928,
    "K": 38.96370668,
    "Mg": 23.9850417,
}

ISOTOPE_ABUNDANCES = {
    "H": [(1.007825032, 0.999885), (2.014101778, 0.000115)],
    "C": [(12.0, 0.9893), (13.003354835, 0.0107)],
    "N": [(14.003074005, 0.99632), (15.000108898, 0.00368)],
    "O": [
        (15.99491462, 0.99757),
        (16.999131757, 0.00038),
        (17.999159613, 0.00205),
    ],
    "S": [
        (31.972071, 0.9499),
        (32.971458, 0.0075),
        (33.967867, 0.0425),
        (35.967081, 0.0001),
    ],
    "P": [(30.973761998, 1.0)],
    "F": [(18.998403163, 1.0)],
    "Cl": [(34.96885268, 0.7576), (36.96590260, 0.2424)],
    "Br": [(78.9183376, 0.5069), (80.916291, 0.4931)],
    "I": [(126.904468, 1.0)],
    "Si": [(27.9769265, 0.92223), (28.9764947, 0.04685), (29.9737701, 0.03092)],
    "B": [(10.012937, 0.199), (11.009305, 0.801)],
    "Se": [
        (73.922476, 0.0089),
        (75.919214, 0.0937),
        (76.919915, 0.0763),
        (77.917309, 0.2377),
        (79.916522, 0.4961),
        (81.916700, 0.0873),
    ],
    "Li": [(6.015123, 0.075), (7.016004, 0.925)],
    "Na": [(22.98976928, 1.0)],
    "K": [(38.96370668, 0.9326), (39.96399848, 0.000117), (40.96182576, 0.0673)],
    "Mg": [(23.9850417, 0.7899), (24.9858369, 0.10), (25.9825929, 0.1101)],
}

ELEMENT_SYMBOLS = {
    "C",
    "N",
    "O",
    "S",
    "P",
    "F",
    "Cl",
    "Br",
    "I",
    "H",
    "Si",
    "B",
    "Se",
    "Li",
    "Na",
    "K",
    "Mg",
}

IMPLICIT_H_ELEMENTS = {"C", "N", "O", "S", "P"}
HETERO_ELECTRON_ATOMS = {"N", "O", "S", "P", "F", "Cl", "Br", "I", "Se", "B"}
ELECTRON_SLOT_ANGLES = [0, 45, 90, 135, 180, 225, 270, 315]

SYMBOL_TEXT_TOOLS = {
    "tool_symbol_plus": {"text": "+", "scale": 1.0, "anchor": True},
    "tool_symbol_minus": {"text": "-", "scale": 1.0, "anchor": True},
    "tool_symbol_radical": {"text": "•", "scale": 1.2, "anchor": True, "rotate": False, "electrons": 1},
    "tool_symbol_lone_pair": {"text": "··", "scale": 1.1, "anchor": True, "rotate": True, "electrons": 2},
    "tool_symbol_radical_cation": {"text": "•+", "scale": 1.1, "anchor": True, "rotate": False},
    "tool_symbol_radical_anion": {"text": "•-", "scale": 1.1, "anchor": True, "rotate": False},
    "tool_symbol_partial_plus": {"text": "δ+", "scale": 1.0, "anchor": True, "rotate": False},
    "tool_symbol_partial_minus": {"text": "δ-", "scale": 1.0, "anchor": True, "rotate": False},
}
