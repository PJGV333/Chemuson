"""Controllers para desacoplar lógica de la ventana principal."""

from .clean2d_controller import Clean2DController
from .document_controller import (
    DocumentController,
    DocumentDiscardContext,
    DocumentTabsContext,
    RecentFilesContext,
)
from .export_controller import ExportController
from .file_controller import FileController, FileWorkflowContext
from .recovery_controller import RecoveryController
from .text_format_controller import TextFormatController
from .template_controller import TemplateController, TemplateControllerContext
from .update_controller import UpdateController, UpdateControllerContext
from .view_controller import ViewController

__all__ = [
    "Clean2DController",
    "DocumentController",
    "DocumentDiscardContext",
    "DocumentTabsContext",
    "ExportController",
    "FileController",
    "FileWorkflowContext",
    "RecentFilesContext",
    "RecoveryController",
    "TextFormatController",
    "TemplateController",
    "TemplateControllerContext",
    "UpdateController",
    "UpdateControllerContext",
    "ViewController",
]
