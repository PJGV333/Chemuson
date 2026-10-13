from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

from PyQt6 import sip
from PyQt6.QtCore import QCoreApplication, QEvent, QEventLoop, QThread, QTimer
from PyQt6.QtWidgets import QApplication, QFileDialog

from chemuson.chemio.rdkit_io import molfile_to_molgraph, smiles_to_molgraph
from chemuson.gui.main_window import ChemusonWindow


SMILES = "C[C@H](F)CO"


def _wait_until_closed(window: ChemusonWindow, timeout_ms: int = 12000) -> bool:
    if not window.isVisible():
        return True
    loop = QEventLoop()
    poll = QTimer()
    poll.setInterval(10)
    poll.timeout.connect(lambda: loop.quit() if not window.isVisible() else None)
    poll.start()
    QTimer.singleShot(timeout_ms, loop.quit)
    loop.exec()
    poll.stop()
    return not window.isVisible()


def _window() -> tuple[QApplication, ChemusonWindow]:
    app = QApplication([])
    window = ChemusonWindow()
    window._maybe_show_onboarding = lambda: None
    window._confirm_discard_changes = lambda _canvas: True
    window._apply_pending_portable_update_on_exit = lambda: True
    window._apply_pending_windows_update_on_exit = lambda: True
    window.show()
    app.processEvents()
    return app, window


def _deferred_delete(app: QApplication, widget) -> None:
    widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()


def main() -> int:
    scenario = sys.argv[1]
    result: dict[str, object] = {"scenario": scenario}
    app, window = _window()
    canvas = window.canvas

    if scenario in {"exported", "modified", "recovered"}:
        canvas._insert_molgraph(smiles_to_molgraph(SMILES))
        app.processEvents()

    if scenario == "exported":
        root = Path(tempfile.mkdtemp(prefix="chemuson-p0-mol-"))
        target = root / "prueba_mol_1"
        QFileDialog.getSaveFileName = lambda *_args, **_kwargs: (
            str(target),
            "Archivo MOL (*.mol)",
        )
        window._on_file_save()
        mol_path = target.with_suffix(".mol")
        assert mol_path.is_file()
        assert not target.exists()
        imported = molfile_to_molgraph(mol_path.read_text(encoding="utf-8"))
        assert len(imported.atoms) == 5 and len(imported.bonds) == 4
        assert window._tab_manager.file_path_for(canvas) == str(mol_path)
        assert canvas.undo_stack.isClean()
        result["exported"] = str(mol_path)

    if scenario == "modified":
        assert not canvas.undo_stack.isClean()

    if scenario == "recovered":
        manager = window._tab_manager.autosave_manager_for(canvas)
        snapshot = manager.core._write_autosave(force=True)
        assert snapshot and Path(snapshot).is_file()
        assert window._open_autosave_document(snapshot, original_path="/tmp/recovered.cmsn")
        canvas = window.canvas
        assert len(canvas.model.atoms) == 5 and len(canvas.model.bonds) == 4
        carbon_labels = [
            item
            for item in canvas.atom_items.values()
            if item.element == "C" and item.label.isVisible()
        ]
        visible_labels = [
            item.label.toPlainText()
            for item in canvas.atom_items.values()
            if item.label.isVisible()
        ]
        assert not carbon_labels, visible_labels
        assert len(visible_labels) == 2, visible_labels
        result["recovered"] = str(snapshot)
        result["atoms"] = len(canvas.model.atoms)
        result["bonds"] = len(canvas.model.bonds)
        result["visible_labels"] = visible_labels

    if scenario == "canvas-delete":
        canvas_ref = canvas
        assert window._close_canvas_tab(canvas_ref)
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
        app.processEvents()
        assert sip.isdeleted(canvas_ref)
        canvas = window.canvas
        result["canvas_deleted"] = True

    if scenario == "async-cleanup":
        class FiniteThread(QThread):
            def run(self) -> None:
                time.sleep(0.05)

        thread = FiniteThread(window)
        thread.start()
        result["thread_started"] = True

    window.close()
    assert window._active_canvas_connected is None
    assert _wait_until_closed(window), "window did not complete close before timeout"
    remaining = [thread for thread in window.findChildren(QThread) if thread.isRunning()]
    assert not remaining, f"running descendant QThreads after close: {len(remaining)}"
    result["closed"] = True
    result["remaining_threads"] = len(remaining)

    if not sip.isdeleted(canvas):
        _deferred_delete(app, canvas)
    _deferred_delete(app, window)
    assert sip.isdeleted(window)
    result["canvas_deleted"] = sip.isdeleted(canvas)
    result["window_deleted"] = True
    print(json.dumps(result, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
