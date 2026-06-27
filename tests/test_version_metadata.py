"""Pruebas de versión única para app y CLI."""

import sys


from chemuson import __version__
from chemuson.__main__ import main
from chemuson.version import get_app_version, get_installed_version


def test_single_version_source_matches_public_api() -> None:
    assert get_app_version() == __version__
    assert isinstance(get_installed_version(), str)


def test_cli_version_flag_prints_version(capsys, monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["chemuson", "--version"])
    main()
    out = capsys.readouterr().out.strip()
    assert out == __version__

