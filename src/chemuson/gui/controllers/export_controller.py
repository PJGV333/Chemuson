from __future__ import annotations

from PyQt6.QtGui import QPainter
from PyQt6.QtPrintSupport import QPrinter
from PyQt6.QtWidgets import QFileDialog, QMessageBox


class ExportController:
    """Delega exportación de imágenes/documentos del canvas activo."""

    def export(self, window, format_name: str) -> None:
        canvas = window.canvas
        if format_name == "png":
            image = canvas._render_scene_image()
            if image:
                filepath, _ = QFileDialog.getSaveFileName(window, "Exportar PNG", "", "Imagen PNG (*.png)")
                if filepath:
                    image.save(filepath, "PNG")
                    window.statusBar().showMessage(f"Exportado: {filepath}")
            return

        if format_name == "svg":
            svg_data = canvas._render_scene_svg()
            filepath, _ = QFileDialog.getSaveFileName(window, "Exportar SVG", "", "Imagen SVG (*.svg)")
            if not filepath:
                return
            if svg_data:
                with open(filepath, "w", encoding="utf-8") as f:
                    f.write(svg_data.decode("utf-8", errors="replace"))
            else:
                from chemuson.chemio.rdkit_io import molgraph_to_svg

                svg_text = molgraph_to_svg(canvas.graph)
                with open(filepath, "w", encoding="utf-8") as f:
                    f.write(svg_text)
                if canvas.state.numbering_enabled and canvas.state.numbering_include_in_export:
                    QMessageBox.information(
                        window,
                        "SVG sin numeración",
                        "Este entorno no permite render SVG desde la escena.\n"
                        "Se exportó SVG químico base sin overlay de numeración.",
                    )
            window.statusBar().showMessage(f"Exportado: {filepath}")
            return

        if format_name == "pdf":
            filepath, _ = QFileDialog.getSaveFileName(window, "Exportar PDF", "", "Documento PDF (*.pdf)")
            if not filepath:
                return
            try:
                printer = QPrinter(QPrinter.PrinterMode.HighResolution)
                printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
                printer.setOutputFileName(filepath)

                bounds = canvas._render_scene_bounds(selected_only=False)
                if bounds is None:
                    QMessageBox.warning(window, "Error", "No hay contenido para exportar.")
                    return

                def render_pdf() -> bool:
                    painter = QPainter(printer)
                    if not painter.isActive():
                        return False
                    canvas.scene.render(
                        painter,
                        printer.pageRect(QPrinter.Unit.Point),
                        bounds,
                    )
                    painter.end()
                    return True

                ok = canvas._with_hidden_render_items(render_pdf)
                if not ok:
                    QMessageBox.warning(window, "Error", "No se pudo exportar PDF.")
                    return
                window.statusBar().showMessage(f"Exportado: {filepath}")
            except Exception as exc:
                QMessageBox.warning(window, "Error", f"No se pudo exportar PDF:\n{exc}")
