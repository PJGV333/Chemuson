"""Caracterización no Qt del entry point y del ciclo de arranque."""

from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from unittest.mock import ANY

import pytest


ROOT = Path(__file__).resolve().parents[1]
MAIN_WINDOW = ROOT / "src" / "chemuson" / "gui" / "main_window.py"


def _run_isolated(code: str, *args: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src")
    return subprocess.run(
        [sys.executable, "-c", code, *args],
        cwd=ROOT,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize("argument", ["--version", "--help"])
def test_non_graphical_cli_actions_do_not_load_gui(argument: str) -> None:
    result = _run_isolated(
        """
import sys
from chemuson.__main__ import main

sys.argv = ["chemuson", sys.argv[1]]
try:
    main()
except SystemExit as exc:
    if exc.code not in (0, None):
        raise

forbidden = sorted(
    name
    for name in sys.modules
    if name == "PyQt6"
    or name.startswith("PyQt6.")
    or name == "chemuson.gui"
    or name.startswith("chemuson.gui.")
)
if forbidden:
    raise AssertionError(forbidden)
""",
        argument,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip()


def test_importing_entry_point_does_not_load_gui() -> None:
    result = _run_isolated(
        """
import sys
import chemuson.__main__

forbidden = sorted(
    name
    for name in sys.modules
    if name == "PyQt6"
    or name.startswith("PyQt6.")
    or name == "chemuson.gui"
    or name.startswith("chemuson.gui.")
)
if forbidden:
    raise AssertionError(forbidden)
"""
    )

    assert result.returncode == 0, result.stderr


def test_default_cli_invocation_delegates_once() -> None:
    result = _run_isolated(
        """
import sys
from types import ModuleType

calls = []
gui = ModuleType("chemuson.gui")
main_window = ModuleType("chemuson.gui.main_window")
main_window.run_app = lambda: calls.append("run_app")
sys.modules["chemuson.gui"] = gui
sys.modules["chemuson.gui.main_window"] = main_window

from chemuson.__main__ import main

sys.argv = ["chemuson"]
main()
if calls != ["run_app"]:
    raise AssertionError(calls)
"""
    )

    assert result.returncode == 0, result.stderr


def _load_historical_run_app(namespace: dict[str, Any]) -> Any:
    tree = ast.parse(MAIN_WINDOW.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "run_app"
    )
    module = ast.Module(body=[function], type_ignores=[])
    ast.fix_missing_locations(module)
    exec(compile(module, MAIN_WINDOW, "exec"), namespace)
    return namespace["run_app"]


class _ExitObserved(Exception):
    def __init__(self, code: int) -> None:
        self.code = code


def test_historical_startup_success_order() -> None:
    events: list[Any] = []

    class FakeApplication:
        _instance: FakeApplication | None = None

        def __init__(self, argv: list[str]) -> None:
            type(self)._instance = self
            events.append(("create_application", argv))

        @classmethod
        def instance(cls) -> FakeApplication | None:
            return cls._instance

        def setApplicationName(self, name: str) -> None:
            events.append(("set_name", name))

        def setApplicationVersion(self, version: str) -> None:
            events.append(("set_version", version))

        def exec(self) -> int:
            events.append("event_loop")
            return 23

    class FakeWindow:
        def __init__(self) -> None:
            events.append("create_window")

        @staticmethod
        def check_autosaves(window: FakeWindow) -> None:
            events.append(("check_autosaves", window))

        def show(self) -> None:
            events.append("show_window")

    namespace = {
        "QApplication": FakeApplication,
        "ChemusonWindow": FakeWindow,
        "QMessageBox": SimpleNamespace(),
        "crash_reporter": SimpleNamespace(
            install=lambda: events.append("install_crash_reporter")
        ),
        "get_app_version": lambda: "1.2.3",
        "sys": SimpleNamespace(
            argv=["chemuson"],
            exit=lambda code: (_ for _ in ()).throw(_ExitObserved(code)),
        ),
    }
    run_app = _load_historical_run_app(namespace)

    with pytest.raises(_ExitObserved) as observed:
        run_app()

    assert observed.value.code == 23
    assert events == [
        "install_crash_reporter",
        ("create_application", ["chemuson"]),
        ("set_name", "Chemuson"),
        ("set_version", "1.2.3"),
        "create_window",
        ("check_autosaves", ANY),
        "show_window",
        "event_loop",
    ]


@pytest.mark.parametrize("application_exists", [False, True])
def test_historical_startup_failure_reporting(application_exists: bool) -> None:
    events: list[Any] = []
    stderr: list[str] = []
    application = object() if application_exists else None

    class FakeApplication:
        def __init__(self, argv: list[str]) -> None:
            raise RuntimeError("startup failed")

        @staticmethod
        def instance() -> object | None:
            return application

    namespace = {
        "QApplication": FakeApplication,
        "ChemusonWindow": object,
        "QMessageBox": SimpleNamespace(
            critical=lambda *args: events.append(("critical", args))
        ),
        "crash_reporter": SimpleNamespace(
            install=lambda: events.append("install_crash_reporter"),
            write_crash_log=lambda exc: events.append(("write_log", str(exc)))
            or Path("/tmp/chemuson-crash.log"),
        ),
        "get_app_version": lambda: "1.2.3",
        "sys": SimpleNamespace(
            argv=["chemuson"],
            stderr=SimpleNamespace(write=stderr.append),
            exit=lambda code: events.append(("exit", code)),
        ),
    }
    run_app = _load_historical_run_app(namespace)

    run_app()

    assert events[:2] == [
        "install_crash_reporter",
        ("write_log", "startup failed"),
    ]
    if application_exists:
        assert events[2][0] == "critical"
        assert "No se pudo iniciar la aplicación." in events[2][1][2]
        assert stderr == []
    else:
        assert events == events[:2]
        assert "".join(stderr) == (
            "No se pudo iniciar la aplicación.\n"
            "Se guardó un reporte en: /tmp/chemuson-crash.log\n"
        )
