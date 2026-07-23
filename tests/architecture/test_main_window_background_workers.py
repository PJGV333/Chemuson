"""Contratos arquitectónicos de los workers internos de la ventana."""

from __future__ import annotations

import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MAIN_WINDOW = ROOT / "src" / "chemuson" / "gui" / "main_window.py"
WORKERS = ROOT / "src" / "chemuson" / "gui" / "background_workers.py"


def _classes(path: Path) -> dict[str, ast.ClassDef]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {
        node.name: node
        for node in tree.body
        if isinstance(node, ast.ClassDef)
    }


def test_workers_have_single_internal_owner() -> None:
    worker_classes = _classes(WORKERS)
    window_classes = _classes(MAIN_WINDOW)

    assert {"_DescriptorWorker", "_NameToStructureWorker"} <= worker_classes.keys()
    assert "_DescriptorWorker" not in window_classes
    assert "_NameToStructureWorker" not in window_classes


def test_main_window_imports_workers_and_retains_thread_orchestration() -> None:
    source = MAIN_WINDOW.read_text(encoding="utf-8")

    assert "from chemuson.gui.background_workers import (" in source
    assert "QThread(self)" in source
    assert "thread.started.connect(worker.run)" in source
    assert "worker.finished.connect(thread.quit)" in source
    assert "_descriptor_jobs[job_id]" in source
    assert "_name2structure_jobs[job_id]" in source


def test_worker_contracts_preserve_signals_and_deferred_backends() -> None:
    source = WORKERS.read_text(encoding="utf-8")

    assert "finished = pyqtSignal(int, dict, str)" in source
    assert "finished = pyqtSignal(int, object, str)" in source
    assert "molecular_descriptors_isolated(self._graph, timeout_s=5.0)" in source
    assert "resolve_name_to_structure(self._query, allow_network=True, timeout_s=8.0)" in source

    tree = ast.parse(source)
    top_level_imports = {
        alias.name
        for node in tree.body
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    top_level_imports.update(
        node.module or ""
        for node in tree.body
        if isinstance(node, ast.ImportFrom)
    )
    assert "chemuson.chemio.rdkit_safe" not in top_level_imports
    assert "chemuson.name2structure" not in top_level_imports
    assert "chemuson.gui.main_window" not in top_level_imports
