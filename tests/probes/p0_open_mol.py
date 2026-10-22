from __future__ import annotations

import json
import sys
from pathlib import Path

from PyQt6 import sip
from PyQt6.QtCore import QCoreApplication, QEvent, QEventLoop, QThread, QTimer
from PyQt6.QtWidgets import QApplication, QMessageBox
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemuson.chemio.rdkit_io import molgraph_to_smiles
from chemuson.gui.main_window import ChemusonWindow


def _canonical(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _wait_for_window_close(window: ChemusonWindow, timeout_ms: int = 12000) -> bool:
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


def main() -> int:
    mol_path = Path(sys.argv[1]).resolve()
    expected_smiles = sys.argv[2]
    app = QApplication([])
    window = ChemusonWindow()
    window._maybe_show_onboarding = lambda: None
    window._confirm_discard_changes = lambda _canvas: True
    window._apply_pending_portable_update_on_exit = lambda: True
    window._apply_pending_windows_update_on_exit = lambda: True
    window._properties_update_timer.stop()
    errors: list[str] = []
    QMessageBox.critical = lambda _parent, _title, message: errors.append(str(message))
    window.show()
    app.processEvents()

    window._open_file_path(str(mol_path))
    app.processEvents()
    assert not errors, errors
    assert window.tabs.count() == 2
    canvas = window.canvas
    assert window._tab_manager.file_path_for(canvas) == str(mol_path)
    assert len(canvas.model.atoms) == 5
    assert len(canvas.model.bonds) == 4
    smiles = molgraph_to_smiles(canvas.model)
    molecule = Chem.MolFromSmiles(smiles)
    expected = Chem.MolFromSmiles(expected_smiles)
    assert molecule is not None and expected is not None
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    assert canonical == _canonical(expected_smiles)
    assert rdMolDescriptors.CalcMolFormula(molecule) == "C3H7FO"
    cip = tuple(
        sorted(
            label
            for _index, label in Chem.FindMolChiralCenters(
                molecule,
                includeUnassigned=False,
                useLegacyImplementation=False,
            )
        )
    )
    expected_cip = tuple(
        sorted(
            label
            for _index, label in Chem.FindMolChiralCenters(
                expected,
                includeUnassigned=False,
                useLegacyImplementation=False,
            )
        )
    )
    assert cip == expected_cip, {"actual_cip": cip, "expected_cip": expected_cip, "smiles": smiles}

    atom_count, bond_count = len(canvas.model.atoms), len(canvas.model.bonds)
    # Keep this probe focused on import; close-lifecycle is covered separately
    # by the shutdown subprocess suite, which intentionally does not disconnect here.
    window._before_canvas_discard(canvas)
    window.close()
    assert _wait_for_window_close(window), "window close did not finish before timeout"
    assert not any(thread.isRunning() for thread in window.findChildren(QThread))
    canvas.deleteLater()
    QCoreApplication.sendPostedEvents(canvas, QEvent.Type.DeferredDelete)
    app.processEvents()
    assert sip.isdeleted(canvas)
    window.deleteLater()
    QCoreApplication.sendPostedEvents(window, QEvent.Type.DeferredDelete)
    app.processEvents()
    assert sip.isdeleted(window)
    print(
        json.dumps(
            {
                "opened": str(mol_path),
                "canonical_isomeric_smiles": canonical,
                "formula": rdMolDescriptors.CalcMolFormula(molecule),
                "atoms": atom_count,
                "bonds": bond_count,
                "cip": cip,
                "closed": True,
                "window_deleted": True,
            },
            sort_keys=True,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
