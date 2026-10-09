"""Fail-closed validation for Chemuson's Linux AppImage Type 2 package."""

from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

APP_ID = "io.github.PJGV333.Chemuson"
_SHA_RE = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$", re.IGNORECASE)


def validate_type2_header(path: Path) -> str:
    """Validate x86_64 ELF plus the AppImage Type 2 marker and return SHA-256."""
    if path.suffix != ".AppImage":
        raise ValueError(f"Linux package must retain the .AppImage suffix: {path.name}")
    try:
        with path.open("rb") as stream:
            header = stream.read(64)
    except OSError as exc:
        raise ValueError(f"Cannot read AppImage package: {path}") from exc
    if len(header) < 20 or header[:4] != b"\x7fELF":
        raise ValueError("AppImage payload is not a complete ELF executable.")
    if header[4] != 2 or header[5] != 1:
        raise ValueError("AppImage must be a little-endian 64-bit ELF executable.")
    if int.from_bytes(header[18:20], "little") != 62:
        raise ValueError("AppImage ELF architecture is not x86_64.")
    if header[8:11] != b"AI\x02":
        raise ValueError("ELF payload lacks the AppImage Type 2 AI\\x02 signature.")

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run(args: list[str], *, cwd: Path, env: dict[str, str], timeout: int) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            args,
            cwd=cwd,
            env=env,
            check=True,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        details = getattr(exc, "stderr", "") or getattr(exc, "stdout", "") or str(exc)
        raise ValueError(f"AppImage validation command failed: {args[0]}: {details}") from exc


def _validate_pyinstaller_resources(executable: Path) -> None:
    try:
        from PyInstaller.archive.readers import CArchiveReader

        names = {
            str(name).replace("\\", "/").casefold()
            for name in CArchiveReader(str(executable)).toc
        }
    except Exception as exc:
        raise ValueError("Could not inspect the bundled PyInstaller CArchive.") from exc

    if not any(name.startswith(("pyqt6/qtcore.", "pyqt6.qtcore")) for name in names):
        raise ValueError("PyInstaller CArchive does not contain the PyQt6 QtCore module.")
    if not any(name.startswith(("pyqt6/qtwidgets.", "pyqt6.qtwidgets")) for name in names):
        raise ValueError("PyInstaller CArchive does not contain the PyQt6 QtWidgets module.")
    if not any(name.startswith(("pyqt6/qtsvg.", "pyqt6.qtsvg")) for name in names):
        raise ValueError("PyInstaller CArchive does not contain the PyQt6 QtSvg module.")
    packaged_svgs = {
        name
        for name in names
        if "chemuson/gui/theme/icons/i-" in name and name.endswith(".svg")
    }
    if len(packaged_svgs) != 69:
        raise ValueError(
            f"PyInstaller CArchive must contain all 69 static SVG icons; found {len(packaged_svgs)}."
        )
    if not any("/plugins/platforms/libqoffscreen.so" in name for name in names):
        raise ValueError("PyInstaller CArchive does not contain Qt's offscreen platform plugin.")
    if not any(name.endswith("chemuson/gui/theme/icons/i-flask.svg") for name in names):
        raise ValueError("PyInstaller CArchive does not contain ChemUSON's packaged SVG resources.")
    expected_templates = {
        "chemuson/chemname/templates/fused/pyrene_cas.mol",
        "chemuson/chemname/templates/fused/pyrene_iupac2004.mol",
        "chemuson/chemname/templates/simple/benzene.mol",
        "chemuson/chemname/templates/special/alpha_d_glucopyranose.mol",
        "chemuson/chemname/templates/special/androstane_core.mol",
        "chemuson/chemname/templates/special/beta_d_fructofuranose.mol",
        "chemuson/chemname/templates/special/beta_d_glucopyranose.mol",
        "chemuson/chemname/templates/special/cholestane_core.mol",
        "chemuson/chemname/templates/special/d_ribose.mol",
    }
    packaged_templates = {
        name for name in names if "/chemuson/chemname/templates/" in f"/{name}"
    }
    if packaged_templates != expected_templates:
        raise ValueError(
            "PyInstaller CArchive must contain exactly the nine ChemName MOL templates; "
            f"found {sorted(packaged_templates)}."
        )


def _validate_appdir(appdir: Path, *, version: str) -> tuple[Path, Path]:
    apprun = appdir / "AppRun"
    desktop = appdir / f"{APP_ID}.desktop"
    icon = appdir / f"{APP_ID}.svg"
    executable = appdir / "usr/bin/Chemuson"
    appstream = appdir / f"usr/share/metainfo/{APP_ID}.appdata.xml"
    installed_desktop = appdir / f"usr/share/applications/{APP_ID}.desktop"
    installed_icon = appdir / f"usr/share/icons/hicolor/scalable/apps/{APP_ID}.svg"

    for required in (apprun, desktop, icon, executable, appstream, installed_desktop, installed_icon):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"AppDir is missing a required nonempty entry: {required.relative_to(appdir)}")
    if not os.access(apprun, os.X_OK) or not os.access(executable, os.X_OK):
        raise ValueError("AppRun and usr/bin/Chemuson must be executable.")
    launcher = apprun.read_text(encoding="utf-8")
    if "usr/bin/Chemuson" not in launcher or "exec " not in launcher:
        raise ValueError("AppRun must exec the packaged Chemuson binary.")

    parser = configparser.ConfigParser(interpolation=None, strict=False)
    parser.optionxform = str
    try:
        parser.read_string(desktop.read_text(encoding="utf-8"))
    except (configparser.Error, UnicodeDecodeError) as exc:
        raise ValueError("AppDir desktop entry is invalid.") from exc
    entry = parser["Desktop Entry"] if parser.has_section("Desktop Entry") else {}
    if entry.get("Type") != "Application" or entry.get("Exec") != "AppRun":
        raise ValueError("AppDir desktop entry must be an Application launching AppRun.")
    if entry.get("Name") != "ChemUSON" or entry.get("Icon") != APP_ID:
        raise ValueError("AppDir desktop entry name/icon does not identify ChemUSON.")
    if entry.get("X-AppImage-Version") != version:
        raise ValueError("AppDir desktop metadata version does not match the canonical version.")

    validator = shutil.which("desktop-file-validate")
    if not validator:
        raise ValueError("desktop-file-validate is required to validate the AppDir desktop entry.")
    for desktop_path in (desktop, installed_desktop):
        _run([validator, str(desktop_path)], cwd=appdir, env=os.environ.copy(), timeout=15)

    try:
        icon_root = ET.parse(icon).getroot()
        installed_icon_root = ET.parse(installed_icon).getroot()
        metadata_root = ET.parse(appstream).getroot()
    except ET.ParseError as exc:
        raise ValueError("AppDir SVG/AppStream resource is not valid XML.") from exc
    if icon_root.tag.rsplit("}", 1)[-1] != "svg" or installed_icon_root.tag.rsplit("}", 1)[-1] != "svg":
        raise ValueError("AppDir application icon must be SVG.")
    metadata_id = metadata_root.findtext("id", default="").strip()
    if metadata_id != APP_ID:
        raise ValueError("AppDir AppStream component ID does not match ChemUSON.")
    releases = metadata_root.findall(".//release")
    if not any(release.attrib.get("version") == version for release in releases):
        raise ValueError("AppDir AppStream release version does not match the canonical version.")
    _validate_pyinstaller_resources(executable)
    return apprun, executable


def _validate_frozen_icons(executable: Path, *, cwd: Path, env: dict[str, str]) -> dict[str, object]:
    validator = Path(__file__).with_name("validate_packaged_icons.py").resolve()
    result = _run(
        [sys.executable, str(validator), "--executable", str(executable), "--timeout", "120"],
        cwd=cwd,
        env=env,
        timeout=150,
    )
    try:
        report = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError("Frozen AppImage executable did not return icon smoke JSON.") from exc
    print("Frozen AppImage SVG icon smoke passed for light/dark themes and HiDPI.")
    return report


def _validate_frozen_rdkit_worker(
    executable: Path, *, cwd: Path, env: dict[str, str]
) -> dict[str, object]:
    validator = Path(__file__).with_name("validate_packaged_rdkit_worker.py").resolve()
    result = _run(
        [sys.executable, str(validator), "--executable", str(executable), "--timeout", "120"],
        cwd=cwd,
        env=env,
        timeout=150,
    )
    try:
        report = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError("Frozen AppImage executable did not return RDKit smoke JSON.") from exc
    if not isinstance(report, dict) or report.get("ok") is not True:
        raise ValueError("Frozen AppImage executable failed the RDKit worker smoke.")
    print("Frozen AppImage RDKit imports, descriptors, SMILES, and 3D worker smoke passed.")
    return report


def _validate_frozen_chemname(
    executable: Path, *, cwd: Path, env: dict[str, str]
) -> dict[str, object]:
    validator = Path(__file__).with_name("validate_packaged_chemname.py").resolve()
    result = _run(
        [sys.executable, str(validator), "--executable", str(executable), "--timeout", "120"],
        cwd=cwd,
        env=env,
        timeout=150,
    )
    try:
        report = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise ValueError("Frozen AppImage executable did not return ChemName smoke JSON.") from exc
    if not isinstance(report, dict) or report.get("ok") is not True:
        raise ValueError("Frozen AppImage executable failed the ChemName smoke.")
    print(
        "Frozen AppImage ChemName templates and molecule names match Python source: "
        f"{len(report.get('template_resources', []))} templates, "
        f"{len(report.get('molecule_results', []))} molecule cases."
    )
    return report


def _appimage_update_information(appimage: Path, *, cwd: Path, env: dict[str, str]) -> str:
    result = _run(
        [str(appimage), "--appimage-updateinformation"],
        cwd=cwd,
        env=env,
        timeout=30,
    )
    return result.stdout.strip()


def _headless_startup_smoke(apprun: Path, *, scratch: Path) -> None:
    home = scratch / "smoke-home"
    config = home / ".config"
    data = home / ".local/share"
    cache = home / ".cache"
    for directory in (home, config, data, cache):
        directory.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(
        {
            "HOME": str(home),
            "XDG_CONFIG_HOME": str(config),
            "XDG_DATA_HOME": str(data),
            "XDG_CACHE_HOME": str(cache),
            "QT_QPA_PLATFORM": "offscreen",
        }
    )
    env.pop("APPIMAGE", None)
    env.pop("APPDIR", None)
    log_path = scratch / "headless-launch.log"
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            [str(apprun)],
            cwd=scratch,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            return_code = process.wait(timeout=12)
        except subprocess.TimeoutExpired:
            return_code = None
        if return_code is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)
        elif return_code == 0:
            raise ValueError("Headless GUI process exited before the controlled startup window.")
        else:
            raise ValueError(f"Headless GUI process exited with code {return_code}.")

    crash_logs = list((config / "chemuson/crash_logs").glob("crash_*.txt"))
    if crash_logs:
        details = crash_logs[0].read_text(encoding="utf-8", errors="replace")[-4000:]
        raise ValueError(f"ChemUSON wrote a crash report during headless startup: {details}")
    print("Headless launch stayed alive for 12s; this is not interactive GUI acceptance.")


def _expected_update_info(owner: str, repo: str, channel: str, appimage_name: str) -> str:
    track = "prerelease" if channel == "beta" else "latest"
    return f"gh-releases-zsync|{owner}|{repo}|{track}|{appimage_name}.zsync"


def validate_appimage(
    *,
    appimage: Path,
    version: str,
    source_sha: str,
    build_type: str,
    channel: str = "",
    tag: str = "",
    owner: str = "PJGV333",
    repo: str = "Chemuson",
) -> dict[str, object]:
    if build_type not in {"preview", "release"}:
        raise ValueError(f"Unsupported AppImage build type: {build_type!r}")
    normalized_sha = str(source_sha or "").strip().lower()
    if not _SHA_RE.fullmatch(normalized_sha):
        raise ValueError("AppImage source provenance requires a full Git SHA.")
    appimage = appimage.resolve()
    package_sha = validate_type2_header(appimage)
    environment = os.environ.copy()

    with tempfile.TemporaryDirectory(prefix="chemuson-appimage-validate-") as temp:
        scratch = Path(temp)
        extraction = _run(
            [str(appimage), "--appimage-extract"],
            cwd=scratch,
            env=environment,
            timeout=300,
        )
        appdir = scratch / "squashfs-root"
        if not appdir.is_dir() or "squashfs-root" not in extraction.stdout:
            raise ValueError("AppImage --appimage-extract did not produce squashfs-root.")
        apprun, executable = _validate_appdir(appdir, version=version)
        icon_report = _validate_frozen_icons(executable, cwd=scratch, env=environment)
        rdkit_report = _validate_frozen_rdkit_worker(executable, cwd=scratch, env=environment)
        chemname_report = _validate_frozen_chemname(executable, cwd=scratch, env=environment)
        version_result = _run(
            [str(apprun), "--version"], cwd=scratch, env=environment, timeout=90
        )
        if version_result.stdout.strip() != version:
            raise ValueError(
                f"AppImage internal version {version_result.stdout.strip()!r} does not match {version!r}."
            )

        embedded_update_info = _appimage_update_information(
            appimage, cwd=scratch, env=environment
        )
        updateinfo_path = Path(f"{appimage}.updateinfo")
        update_json_path = Path(f"{appimage}.update.json")
        zsync_path = Path(f"{appimage}.zsync")
        if build_type == "preview":
            if embedded_update_info:
                raise ValueError("Preview AppImage must not embed public updater information.")
            for sidecar in (updateinfo_path, update_json_path, zsync_path):
                if sidecar.exists():
                    raise ValueError(f"Preview AppImage has forbidden public updater sidecar: {sidecar.name}")
        else:
            if channel not in {"beta", "stable"} or not tag:
                raise ValueError("Release AppImage requires a validated channel and tag.")
            expected = _expected_update_info(owner, repo, channel, appimage.name)
            if embedded_update_info != expected:
                raise ValueError("Embedded AppImage update information differs from the existing channel contract.")
            if not updateinfo_path.is_file() or updateinfo_path.read_text(encoding="utf-8").strip() != expected:
                raise ValueError("AppImage .updateinfo sidecar differs from the embedded update information.")
            if not update_json_path.is_file() or update_json_path.stat().st_size == 0:
                raise ValueError("Release AppImage .update.json sidecar is missing or empty.")
            update = json.loads(update_json_path.read_text(encoding="utf-8"))
            expected_fields = {
                "version": version,
                "channel": channel,
                "tag": tag,
                "source_sha": normalized_sha,
                "appimage_update_information": expected,
            }
            for key, value in expected_fields.items():
                if update.get(key) != value:
                    raise ValueError(f"AppImage updater metadata {key!r} does not match its release identity.")
            if not zsync_path.is_file() or zsync_path.stat().st_size == 0:
                raise ValueError("Release AppImage .zsync update payload is missing or empty.")
            zsync_header = zsync_path.read_text(encoding="utf-8", errors="replace").splitlines()
            zsync_url = next((line[5:] for line in zsync_header if line.startswith("URL: ")), "")
            if not zsync_url.endswith(f"/{appimage.name}"):
                raise ValueError("AppImage .zsync URL does not target the contractual release asset.")

        _headless_startup_smoke(apprun, scratch=scratch)

    return {
        "appimage": str(appimage),
        "version": version,
        "source_sha": normalized_sha,
        "sha256": package_sha,
        "build_type": build_type,
        "embedded_update_information": embedded_update_info,
        "icon_smoke": icon_report,
        "rdkit_worker_smoke": rdkit_report,
        "chemname_smoke": chemname_report,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--appimage", type=Path)
    parser.add_argument("--header-only", action="store_true")
    parser.add_argument("--version")
    parser.add_argument("--source-sha")
    parser.add_argument("--build-type", choices=("preview", "release"))
    parser.add_argument("--channel", default="")
    parser.add_argument("--tag", default="")
    parser.add_argument("--owner", default="PJGV333")
    parser.add_argument("--repo", default="Chemuson")
    args = parser.parse_args()
    if args.appimage is None:
        parser.error("--appimage is required")
    if args.header_only:
        digest = validate_type2_header(args.appimage.resolve())
        print(f"Validated AppImage Type 2 ELF header SHA-256: {digest}")
        return
    if not args.version or not args.source_sha or not args.build_type:
        parser.error("--version, --source-sha, and --build-type are required for full validation")
    record = validate_appimage(
        appimage=args.appimage,
        version=args.version,
        source_sha=args.source_sha,
        build_type=args.build_type,
        channel=args.channel,
        tag=args.tag,
        owner=args.owner,
        repo=args.repo,
    )
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"AppImage validation failed: {exc}") from exc
