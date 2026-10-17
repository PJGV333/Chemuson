from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtWidgets import QApplication

from chemuson.gui.style import CHEMDOODLE_LIKE
from tools.visual.depiction_lab_app.named_presets import EDITORIAL_AB3, EDITORIAL_AB3_BASED_ON
from tools.visual.depiction_lab_app.profiles import load_profile, profile_directory, safe_profile_filename
from tools.visual.depiction_lab_app.scenes import CASES, graph_snapshot
from tools.visual.depiction_lab_app.ui import DepictionLabWindow

app = QApplication([])
window = DepictionLabWindow()
window.show()
app.processEvents()

current_graph = graph_snapshot(window.current_panel.adapter.canvas.model)
current_wedge = next(
    item
    for item in window.current_panel.adapter.canvas.bond_items.values()
    if item.style.value == "wedge"
)
working_wedge = next(
    item
    for item in window.working_panel.adapter.canvas.bond_items.values()
    if item.style.value == "wedge"
)
current_path = current_wedge.path()
working_before = working_wedge.path()
control = window.parameter_controls["solid_wedge_width_ratio"]
control.spin.setValue(0.34)
app.processEvents()
assert window.working_parameters.solid_wedge_width_ratio == 0.34
assert current_wedge.path() == current_path
assert working_wedge.path() != working_before
assert graph_snapshot(window.current_panel.adapter.canvas.model) == current_graph

window.zoom_combo.setCurrentIndex(5)
window.theme_combo.setCurrentIndex(1)
window.angle_combo.setCurrentIndex(3)
window.diagnostics_toggle.setChecked(True)
app.processEvents()
assert window.zoom == 3.0 and window.theme == "dark" and window.angle == 60.0
assert window.current_panel.adapter.canvas.transform() == window.working_panel.adapter.canvas.transform()
assert window.current_panel.adapter.diagnostic_items
assert window.working_panel.adapter.diagnostic_items
window.hide_controls_button.setChecked(True)
assert not window.controls.isVisible()
window.hide_controls_button.setChecked(False)

for index, case in enumerate(CASES):
    window.case_combo.setCurrentIndex(index)
    app.processEvents()
    assert window.case_id == case.case_id
    if case.case_id == "B4":
        assert "GEOMETRIC FIXTURE" in window.case_combo.currentText()
        assert "non-editorial" in window.status.text()
        assert "excluded from preset selection / style calibration" in window.status.text()
    window._assert_panels_match()
    assert graph_snapshot(window.current_panel.adapter.canvas.model) == graph_snapshot(
        window.working_panel.adapter.canvas.model
    )

window.editorial_preset_button.click()
app.processEvents()
assert window.working_parameters == EDITORIAL_AB3.parameters
assert window.working_based_on == EDITORIAL_AB3_BASED_ON
assert window.current_panel.adapter.canvas.drawing_style == CHEMDOODLE_LIKE
assert window.working_panel.adapter.canvas.drawing_style == EDITORIAL_AB3.parameters.to_drawing_style()
assert "B4 está excluido" in window.status.text()

window.profile_name.setText("visual-lab-probe")
window._save_profile()
profile_path = profile_directory() / safe_profile_filename("visual-lab-probe")
name, restored = load_profile(profile_path)
assert name == "visual-lab-probe"
assert restored == window.working_parameters
assert json.loads(profile_path.read_text(encoding="utf-8"))["based_on"] == EDITORIAL_AB3_BASED_ON
assert window.status.text().startswith("Perfil guardado")

window.close()
app.processEvents()
assert not window.isVisible()
print(json.dumps({"cases": len(CASES), "working": True, "closed": True}))
