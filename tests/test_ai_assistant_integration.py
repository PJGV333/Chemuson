from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


def _requirements(path: Path) -> set[str]:
    return {
        line.strip().split("[", 1)[0].lower().replace("_", "-")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }


def test_requirements_matches_pyproject_runtime_dependencies():
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    dependency_block = re.search(r"dependencies\s*=\s*\[(.*?)\]", pyproject, re.DOTALL)
    assert dependency_block is not None
    declared = {
        value.lower().replace("_", "-")
        for value in re.findall(r'"([^"]+)"', dependency_block.group(1))
    }
    assert _requirements(ROOT / "requirements.txt") == declared
    assert "pillow" in declared


def test_development_requirements_are_separate_and_readme_documents_fresh_checkout():
    assert _requirements(ROOT / "requirements-dev.txt") == {"pytest", "ruff", "pyyaml"}
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "python -m venv .venv" in readme
    assert "python -m pip install -e ." in readme
    assert "python -m pip install -r requirements-dev.txt" in readme
    assert "activate.fish" in readme
    assert "Open Babel" in readme
