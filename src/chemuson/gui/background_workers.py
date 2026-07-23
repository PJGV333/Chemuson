"""Workers Qt internos para operaciones bloqueantes iniciadas por la GUI."""

from PyQt6.QtCore import QObject, pyqtSignal, pyqtSlot


class _DescriptorWorker(QObject):
    """Worker aislado para descriptores RDKit del dock químico."""

    finished = pyqtSignal(int, dict, str)

    def __init__(self, job_id: int, graph) -> None:
        super().__init__()
        self._job_id = int(job_id)
        self._graph = graph

    @pyqtSlot()
    def run(self) -> None:
        try:
            from chemuson.chemio.rdkit_safe import molecular_descriptors_isolated

            descriptors, error = molecular_descriptors_isolated(self._graph, timeout_s=5.0)
            if error:
                self.finished.emit(self._job_id, {}, str(error))
                return
            self.finished.emit(self._job_id, dict(descriptors or {}), "")
        except Exception as exc:
            self.finished.emit(self._job_id, {}, str(exc))


class _NameToStructureWorker(QObject):
    """Worker para resolver Name->Structure sin bloquear la UI."""

    finished = pyqtSignal(int, object, str)

    def __init__(self, job_id: int, query: str) -> None:
        super().__init__()
        self._job_id = int(job_id)
        self._query = str(query or "").strip()

    @pyqtSlot()
    def run(self) -> None:
        try:
            from chemuson.name2structure import resolve_name_to_structure

            result = resolve_name_to_structure(self._query, allow_network=True, timeout_s=8.0)
            self.finished.emit(self._job_id, result, "")
        except Exception as exc:
            self.finished.emit(self._job_id, None, str(exc))
