from __future__ import annotations

from PyQt6.QtCore import QPointF

from chemuson.core.model import MolGraph
from chemuson.gui.controllers.template_controller import (
    TemplateController,
    TemplateControllerContext,
)


def _simple_graph() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    return graph


class _FakeLibrary:
    def __init__(self, graph: MolGraph | None = None, error: Exception | None = None) -> None:
        self.graph = graph or _simple_graph()
        self.error = error

    def get_template(self, template_id: str) -> dict:
        if self.error is not None:
            raise self.error
        return {"id": template_id, "name": "  Plantilla prueba  "}

    def graph_from_template(self, _template_id: str) -> MolGraph:
        if self.error is not None:
            raise self.error
        return self.graph


class _FakeCanvas:
    def __init__(self) -> None:
        self._last_scene_pos = QPointF(10.0, 20.0)
        self.inserted: list[tuple[MolGraph, QPointF]] = []
        self.cancelled = False
        self._fallback_point = QPointF(30.0, 40.0)

    def _insert_molgraph_at(self, graph: MolGraph, target: QPointF) -> None:
        self.inserted.append((graph, target))

    def cancel_template_insert_mode(self) -> None:
        self.cancelled = True

    def mapToScene(self, _point) -> QPointF:
        return self._fallback_point

    def viewport(self):
        return self

    def rect(self):
        return self

    def center(self):
        return object()


def _context(library: _FakeLibrary, canvas: _FakeCanvas):
    statuses: list[str] = []
    inserted: list[tuple[str, MolGraph]] = []
    return (
        TemplateControllerContext(
            parent=None,
            canvas=canvas,
            template_library=library,  # type: ignore[arg-type]
            show_status=statuses.append,
            refresh_template_views=lambda: None,
            insert_template=lambda label, graph: inserted.append((label, graph)),
        ),
        statuses,
        inserted,
    )


def test_start_template_insert_by_id_enters_click_placement() -> None:
    graph = _simple_graph()
    canvas = _FakeCanvas()
    context, statuses, inserted = _context(_FakeLibrary(graph), canvas)

    TemplateController().start_template_insert_by_id(context, "tpl_1")

    assert inserted == [("Plantilla prueba", graph)]
    assert statuses == []
    assert canvas.inserted == []
    assert not canvas.cancelled


def test_start_template_insert_by_id_places_immediately_at_last_position() -> None:
    graph = _simple_graph()
    canvas = _FakeCanvas()
    context, statuses, inserted = _context(_FakeLibrary(graph), canvas)

    TemplateController().start_template_insert_by_id(context, "tpl_1", place_now=True)

    assert inserted == []
    assert canvas.inserted == [(graph, QPointF(10.0, 20.0))]
    assert canvas.cancelled
    assert statuses == ["Plantilla 'Plantilla prueba' insertada"]


def test_start_template_insert_by_id_places_immediately_at_fallback_position() -> None:
    graph = _simple_graph()
    canvas = _FakeCanvas()
    canvas._last_scene_pos = None
    context, _statuses, _inserted = _context(_FakeLibrary(graph), canvas)

    TemplateController().start_template_insert_by_id(context, "tpl_1", place_now=True)

    assert canvas.inserted == [(graph, QPointF(30.0, 40.0))]


def test_start_template_insert_by_id_reports_load_error(monkeypatch) -> None:
    canvas = _FakeCanvas()
    context, statuses, inserted = _context(_FakeLibrary(error=ValueError("rota")), canvas)
    errors: list[tuple[object, str, str]] = []
    monkeypatch.setattr(
        "chemuson.gui.controllers.template_controller.QMessageBox.critical",
        lambda parent, title, message: errors.append((parent, title, message)),
    )

    TemplateController().start_template_insert_by_id(context, "tpl_rota")

    assert inserted == []
    assert canvas.inserted == []
    assert statuses == []
    assert errors == [(None, "Error", "No se pudo cargar la plantilla:\nrota")]
