from __future__ import annotations

import threading
import time

import pytest
from PyQt6 import sip
from PyQt6.QtCore import QCoreApplication, QEvent
from PyQt6.QtWidgets import QApplication

from chemuson.core.model import MolGraph
from chemuson.gui.dialogs import MolecularAssistantDialog
from chemuson.gui.main_window import ChemusonWindow
from chemuson.molecular_assistant import MolecularAssistantResult, MolecularAssistantStatus


@pytest.fixture(scope="module", autouse=True)
def _qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _isolated_config_home(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))


def _process_qt_deletes() -> None:
    QCoreApplication.processEvents()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    QCoreApplication.processEvents()


def _wait_for(predicate, timeout_s: float = 5.0) -> bool:
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        QCoreApplication.processEvents()
        if predicate():
            return True
        time.sleep(0.01)
    return False


def _success() -> MolecularAssistantResult:
    graph = MolGraph()
    graph.add_atom("C", 0.0, 0.0)
    return MolecularAssistantResult(
        status=MolecularAssistantStatus.SUCCESS,
        provider_id="fake-provider",
        model_id="fake-model",
        proposed_smiles="C",
        graph=graph,
        validation_passed=True,
    )


def _dialog(window: ChemusonWindow) -> MolecularAssistantDialog:
    dialogs = window.findChildren(MolecularAssistantDialog)
    assert dialogs
    return dialogs[-1]


def _open_assistant(window: ChemusonWindow) -> MolecularAssistantDialog:
    window._on_ai_molecular_assistant()
    QApplication.processEvents()
    return _dialog(window)


def _start_job(
    window: ChemusonWindow,
    dialog: MolecularAssistantDialog,
    *,
    generator,
    transform_context=None,
) -> int:
    controller = window._molecular_assistant_controller
    controller._generator = generator
    controller._identity_verifier = None
    window._start_molecular_assistant_job(
        dialog,
        "draw a carbon atom",
        "llama-cpp",
        "http://127.0.0.1:8080/v1",
        "fake-model",
        "",
        False,
        10,
        64,
        transform_context=transform_context,
        identity_enabled=False,
    )
    assert dialog.job_id is not None
    return dialog.job_id


def _assert_empty_registries(window: ChemusonWindow) -> None:
    assert window._molecular_assistant_dialogs == {}
    assert window._molecular_assistant_results == {}
    assert window._molecular_assistant_identity_results == {}
    assert window._molecular_assistant_transform_jobs == {}


def _close_window(window: ChemusonWindow) -> None:
    if not sip.isdeleted(window):
        window.close()
        _process_qt_deletes()


def _close_and_wait(window: ChemusonWindow) -> None:
    accepted = window.close()
    if not accepted:
        assert window._shutdown_started
        assert _wait_for(lambda: window._async_shutdown_complete)
    _process_qt_deletes()
    assert not window.isVisible()


def test_unstarted_dialog_close_and_window_close_leave_no_registry_entries():
    window = ChemusonWindow()
    try:
        dialog = _open_assistant(window)
        dialog.close()
        _process_qt_deletes()

        assert sip.isdeleted(dialog)
        _assert_empty_registries(window)
        assert window.close()
    finally:
        _close_window(window)


def test_pending_dialog_close_abandons_job_and_ignores_late_worker_result():
    window = ChemusonWindow()
    started = threading.Event()
    release = threading.Event()

    def delayed_generator(*_args):
        started.set()
        release.wait(timeout=4.0)
        return _success()

    try:
        dialog = _open_assistant(window)
        job_id = _start_job(window, dialog, generator=delayed_generator)
        assert started.wait(timeout=2.0)

        dialog.close()
        _process_qt_deletes()
        assert sip.isdeleted(dialog)
        _assert_empty_registries(window)
        assert job_id in window._molecular_assistant_controller.active_jobs()
        assert job_id in window._molecular_assistant_controller._abandoned_jobs

        release.set()
        assert _wait_for(lambda: not window._molecular_assistant_controller.active_jobs())
        assert job_id not in window._molecular_assistant_controller._abandoned_jobs
        assert window.close()
    finally:
        release.set()
        _close_window(window)


def test_deleted_qobject_is_removed_before_shutdown_without_close_call():
    """Exercise the real stale-wrapper state: deleteLater without QDialog.finished."""
    window = ChemusonWindow()
    started = threading.Event()
    release = threading.Event()

    def delayed_generator(*_args):
        started.set()
        release.wait(timeout=4.0)
        return _success()

    try:
        dialog = _open_assistant(window)
        job_id = _start_job(window, dialog, generator=delayed_generator)
        assert started.wait(timeout=2.0)

        dialog.deleteLater()
        _process_qt_deletes()
        assert sip.isdeleted(dialog)
        assert job_id not in window._molecular_assistant_dialogs
        assert job_id in window._molecular_assistant_controller._abandoned_jobs
        assert job_id not in window._molecular_assistant_results
        assert job_id not in window._molecular_assistant_identity_results
        assert job_id not in window._molecular_assistant_transform_jobs

        # Before the fix the registry retained the dead wrapper and this close
        # path called dialog.close(), raising PyQt's deleted-C++-object RuntimeError.
        assert not window.close()  # shutdown waits for the deliberately blocked worker
        assert window._shutdown_started
        _assert_empty_registries(window)
        release.set()
        assert _wait_for(lambda: window._async_shutdown_complete)
        _process_qt_deletes()
        assert not window.isVisible()
    finally:
        release.set()
        _close_window(window)


def test_completed_preview_close_removes_result_identity_and_dialog_state():
    window = ChemusonWindow()
    try:
        dialog = _open_assistant(window)
        job_id = _start_job(window, dialog, generator=lambda *_args: _success())
        assert _wait_for(lambda: job_id in window._molecular_assistant_results)
        assert dialog.insert_button.isVisible()
        window._molecular_assistant_identity_results[job_id] = object()

        dialog.close()
        _process_qt_deletes()

        assert sip.isdeleted(dialog)
        _assert_empty_registries(window)
        assert window.close()
    finally:
        _close_window(window)


def test_accept_and_delete_on_close_clean_job_before_application_close():
    window = ChemusonWindow()
    try:
        dialog = _open_assistant(window)
        job_id = _start_job(window, dialog, generator=lambda *_args: _success())
        assert _wait_for(lambda: job_id in window._molecular_assistant_results)

        dialog.insert_button.click()
        _process_qt_deletes()

        assert sip.isdeleted(dialog)
        assert job_id not in window._molecular_assistant_dialogs
        assert job_id not in window._molecular_assistant_results
        assert job_id not in window._molecular_assistant_identity_results
        assert job_id not in window._molecular_assistant_transform_jobs
        assert len(window.canvas.model.atoms) == 1
        window.canvas.undo_stack.setClean()
        assert window.close()
    finally:
        _close_window(window)


def test_repeated_unstarted_open_close_does_not_accumulate_jobs():
    window = ChemusonWindow()
    try:
        for _ in range(5):
            dialog = _open_assistant(window)
            dialog.close()
            _process_qt_deletes()
            assert sip.isdeleted(dialog)
            _assert_empty_registries(window)
        assert window.close()
    finally:
        _close_window(window)


def test_transform_dialog_close_removes_transform_context():
    window = ChemusonWindow()
    started = threading.Event()
    release = threading.Event()

    def delayed_generator(*_args):
        started.set()
        release.wait(timeout=4.0)
        return _success()

    try:
        window._properties_update_timer.stop()
        window.canvas._insert_molgraph(_success().graph, select_inserted=True)
        window.canvas.undo_stack.setClean()
        transform_context = window._capture_selected_molecule_for_ai_transform()
        assert transform_context is not None
        window._open_molecular_assistant_dialog(transform_context)
        dialog = _dialog(window)
        job_id = _start_job(
            window,
            dialog,
            generator=delayed_generator,
            transform_context=transform_context,
        )
        assert started.wait(timeout=2.0)
        assert job_id in window._molecular_assistant_transform_jobs

        dialog.close()
        _process_qt_deletes()
        assert sip.isdeleted(dialog)
        _assert_empty_registries(window)
        release.set()
        assert _wait_for(lambda: not window._molecular_assistant_controller.active_jobs())
        _close_and_wait(window)
    finally:
        release.set()
        _close_window(window)


def test_failed_job_can_retry_with_a_new_id_without_stale_registry_state():
    window = ChemusonWindow()
    failure = MolecularAssistantResult(
        status=MolecularAssistantStatus.PROVIDER_ERROR,
        provider_id="fake-provider",
        model_id="fake-model",
        reason_code="network_error",
    )
    results = [failure, _success()]

    def retry_generator(*_args):
        return results.pop(0)

    try:
        dialog = _open_assistant(window)
        first_job_id = _start_job(window, dialog, generator=retry_generator)
        assert _wait_for(
            lambda: (
                first_job_id not in window._molecular_assistant_dialogs
                and not window._molecular_assistant_controller.active_jobs()
            )
        )
        assert dialog.status_label.text()
        _assert_empty_registries(window)

        second_job_id = _start_job(window, dialog, generator=retry_generator)
        assert second_job_id != first_job_id
        assert _wait_for(lambda: second_job_id in window._molecular_assistant_results)
        assert window._molecular_assistant_dialogs[second_job_id][0] is dialog
        dialog.close()
        _process_qt_deletes()
        assert sip.isdeleted(dialog)
        _assert_empty_registries(window)
        assert window.close()
    finally:
        _close_window(window)


def test_window_shutdown_with_live_dialog_and_worker_waits_safely():
    window = ChemusonWindow()
    started = threading.Event()
    release = threading.Event()

    def delayed_generator(*_args):
        started.set()
        release.wait(timeout=4.0)
        return _success()

    try:
        dialog = _open_assistant(window)
        job_id = _start_job(window, dialog, generator=delayed_generator)
        assert started.wait(timeout=2.0)

        assert not window.close()
        assert window._shutdown_started
        _assert_empty_registries(window)
        release.set()
        assert _wait_for(lambda: window._async_shutdown_complete)
        _process_qt_deletes()
        assert sip.isdeleted(dialog)
        assert not window._molecular_assistant_controller.active_jobs()
        assert not window.isVisible()
    finally:
        release.set()
        _close_window(window)
