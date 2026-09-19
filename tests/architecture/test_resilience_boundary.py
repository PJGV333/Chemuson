from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_resilience_modules_are_canonical_and_shims_are_import_only() -> None:
    for name in ("autosave", "crash_reporter"):
        canonical = ROOT / "src" / "chemuson" / "resilience" / f"{name}.py"
        shim = ROOT / "src" / "chemuson" / "utils" / f"{name}.py"
        assert canonical.exists()
        tree = ast.parse(shim.read_text(encoding="utf-8"))
        assert not any(isinstance(node, (ast.FunctionDef, ast.ClassDef)) for node in tree.body)
        assert "chemuson.resilience" in shim.read_text(encoding="utf-8")


def test_resilience_has_no_chemuson_gui_imports() -> None:
    resilience = ROOT / "src" / "chemuson" / "resilience"
    for path in resilience.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            assert all(not name.startswith("chemuson.gui") for name in names)
