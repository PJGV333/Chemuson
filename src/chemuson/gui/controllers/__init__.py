"""Controllers para desacoplar lógica de la ventana principal."""

from .clean2d_controller import Clean2DController
from .document_controller import DocumentController, DocumentDiscardContext, RecentFilesContext
from .export_controller import ExportController
from .recovery_controller import RecoveryController
from .text_format_controller import TextFormatController
from .template_controller import TemplateController, TemplateControllerContext

__all__ = [
    "Clean2DController",
    "DocumentController",
    "DocumentDiscardContext",
    "ExportController",
    "RecentFilesContext",
    "RecoveryController",
    "TextFormatController",
    "TemplateController",
    "TemplateControllerContext",
]
