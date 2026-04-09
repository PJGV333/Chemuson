"""Harness de acceptance tests para el motor de nomenclatura ChemName."""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import traceback
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from chemuson.chemio.rdkit_io import molfile_to_molgraph, smiles_to_molgraph
from chemuson.chemio.rdkit_safe import (
    is_rdkit_worker_available,
    smiles_to_molgraph_isolated,
)
from chemuson.chemname.engine import iupac_name
from chemuson.chemname.errors import ChemNameNotSupported
from chemuson.chemname.options import NameOptions

DEFAULT_CASES_PATH = Path("tests/data/chemname_acceptance_cases.yml")
DEFAULT_REPORT_PATH = Path("artifacts/chemname_acceptance_report.json")


def load_cases(path: str | Path) -> list[dict[str, Any]]:
    """Carga casos desde YAML/JSON.

    Nota: JSON es un subconjunto de YAML, por lo que también se admite un
    archivo `.yml` con sintaxis JSON.
    """
    cases_path = Path(path)
    text = cases_path.read_text(encoding="utf-8")
    payload: Any
    try:
        payload = json.loads(text)
    except Exception:
        try:
            import yaml  # type: ignore

            payload = yaml.safe_load(text)
        except Exception as exc:
            raise RuntimeError(
                f"No se pudo parsear {cases_path}. Instala PyYAML o usa sintaxis JSON válida."
            ) from exc
    if isinstance(payload, list):
        raw_cases = payload
    elif isinstance(payload, dict):
        raw_cases = payload.get("cases", [])
    else:
        raise ValueError("Formato inválido: se esperaba lista o dict con 'cases'.")

    cases: list[dict[str, Any]] = []
    for idx, raw in enumerate(raw_cases):
        if not isinstance(raw, dict):
            continue
        case = dict(raw)
        case_id = str(case.get("id", "")).strip() or f"case_{idx + 1:03d}"
        case["id"] = case_id
        case["title"] = str(case.get("title", case_id))
        case["input_type"] = str(case.get("input_type", "smiles")).strip().lower()
        case["input"] = str(case.get("input", ""))
        case["expect"] = case.get("expect", r"^(?!N/D$).+")
        case["skip_if_no_rdkit"] = bool(case.get("skip_if_no_rdkit", False))
        case["options"] = dict(case.get("options", {}))
        if "max_ms" in case and case["max_ms"] is not None:
            try:
                case["max_ms"] = float(case["max_ms"])
            except Exception:
                case["max_ms"] = None
        else:
            case["max_ms"] = None
        cases.append(case)
    return cases


def _resolve_file_path(cases_path: Path, raw_path: str) -> Path:
    value = str(raw_path).strip()
    candidate = Path(value)
    if candidate.is_absolute():
        return candidate
    from_cases = (cases_path.parent / candidate).resolve()
    if from_cases.exists():
        return from_cases
    from_cwd = (Path.cwd() / candidate).resolve()
    if from_cwd.exists():
        return from_cwd
    return from_cases


def _looks_like_molblock(text: str) -> bool:
    probe = str(text or "")
    return ("V2000" in probe) or ("V3000" in probe) or ("M  END" in probe)


def _build_options(raw_options: dict[str, Any]) -> NameOptions:
    opts = NameOptions()
    valid = set(getattr(opts, "__dataclass_fields__", {}).keys())
    for key, value in dict(raw_options or {}).items():
        if key not in valid:
            raise ValueError(f"Opción desconocida: {key}")
        setattr(opts, key, value)
    return opts


def _validate_expectation(
    expect: Any,
    *,
    name: str | None,
    exception: Exception | None,
) -> tuple[bool, str]:
    if expect is None:
        expect = r"^(?!N/D$).+"
    expect_text = str(expect).strip()
    resolved_name = str(name or "")
    if expect_text.upper() == "ND":
        if resolved_name == "N/D":
            return True, "expected_nd"
        if isinstance(exception, ChemNameNotSupported):
            return True, "expected_not_supported"
        return False, f"expected ND/NotSupported, got '{resolved_name or 'empty'}'"
    if expect_text.upper() == "NON_ND":
        if resolved_name != "N/D" and resolved_name:
            return True, "expected_non_nd"
        return False, "expected non-N/D"
    try:
        pattern = re.compile(expect_text, re.IGNORECASE)
    except Exception as exc:
        return False, f"regex inválida: {exc}"
    if not resolved_name:
        return False, "empty_name"
    if not pattern.search(resolved_name):
        normalized_name = resolved_name.replace("/", "")
        if normalized_name != resolved_name and pattern.search(normalized_name):
            return True, "regex_match_normalized"
        return False, f"regex '{expect_text}' no coincide con '{resolved_name}'"
    return True, "regex_match"


def _graph_from_smiles(
    smiles: str,
    *,
    opts: NameOptions,
    rdkit_available: bool,
    skip_if_no_rdkit: bool,
) -> tuple[Any | None, str | None, str | None]:
    if opts.rdkit_isolated:
        if not rdkit_available:
            if skip_if_no_rdkit:
                return None, "skip", "rdkit_worker_unavailable"
            return None, "worker_error", "rdkit_worker_unavailable"
        graph, err = smiles_to_molgraph_isolated(smiles)
        if graph is not None:
            return graph, None, None
        return None, "worker_error", str(err or "worker_error")
    try:
        return smiles_to_molgraph(smiles), None, None
    except Exception as exc:
        message = str(exc)
        if ("RDKit no disponible" in message or "No module named" in message) and skip_if_no_rdkit:
            return None, "skip", "rdkit_unavailable"
        return None, "build_error", message


def build_graph_for_case(
    case: dict[str, Any],
    cases_path: str | Path,
    *,
    rdkit_available: bool,
) -> tuple[Any | None, str | None, str | None]:
    """Construye `MolGraph` a partir del caso.

    Returns:
        `(graph, status, reason)`; cuando `status` es `None`, `graph` está listo.
    """
    case_path = Path(cases_path)
    input_type = str(case.get("input_type", "")).strip().lower()
    value = str(case.get("input", "") or "")
    skip_if_no_rdkit = bool(case.get("skip_if_no_rdkit", False))
    opts = _build_options(dict(case.get("options", {})))

    if input_type == "smiles":
        return _graph_from_smiles(
            value,
            opts=opts,
            rdkit_available=rdkit_available,
            skip_if_no_rdkit=skip_if_no_rdkit,
        )

    if input_type == "molblock":
        try:
            return molfile_to_molgraph(value), None, None
        except Exception as exc:
            return None, "build_error", str(exc)

    if input_type == "file":
        filepath = _resolve_file_path(case_path, value)
        if not filepath.exists():
            return None, "build_error", f"file_not_found:{filepath}"
        text = filepath.read_text(encoding="utf-8", errors="replace")
        if _looks_like_molblock(text) or filepath.suffix.lower() in {".mol", ".sdf"}:
            try:
                return molfile_to_molgraph(text), None, None
            except Exception as exc:
                return None, "build_error", str(exc)
        return _graph_from_smiles(
            text.strip(),
            opts=opts,
            rdkit_available=rdkit_available,
            skip_if_no_rdkit=skip_if_no_rdkit,
        )

    return None, "build_error", f"input_type_unsupported:{input_type}"


def evaluate_case(
    case: dict[str, Any],
    cases_path: str | Path,
    *,
    rdkit_available: bool | None = None,
) -> dict[str, Any]:
    """Ejecuta un caso de aceptación y devuelve resultado serializable."""
    if rdkit_available is None:
        rdkit_available = is_rdkit_worker_available()

    started = time.perf_counter()
    build_ms = 0.0
    name_ms = 0.0
    name: str | None = None
    exception: Exception | None = None

    result: dict[str, Any] = {
        "id": str(case.get("id", "")),
        "title": str(case.get("title", "")),
        "status": "error",
        "reason": "",
        "name": "",
        "expect": case.get("expect"),
        "build_ms": 0.0,
        "name_ms": 0.0,
        "duration_ms": 0.0,
        "max_ms": case.get("max_ms"),
        "rdkit_required": bool(case.get("skip_if_no_rdkit", False)),
    }

    try:
        opts = _build_options(dict(case.get("options", {})))
    except Exception as exc:
        result["status"] = "error"
        result["reason"] = f"invalid_options:{exc}"
        return result

    build_start = time.perf_counter()
    graph, build_status, build_reason = build_graph_for_case(
        case,
        cases_path,
        rdkit_available=bool(rdkit_available),
    )
    build_ms = (time.perf_counter() - build_start) * 1000.0

    if build_status == "skip":
        result["status"] = "skip"
        result["reason"] = str(build_reason or "skipped")
    elif build_status == "worker_error":
        result["worker_error"] = str(build_reason or "worker_error")
        if str(case.get("expect", "")).strip().upper() == "ND":
            name = "N/D"
            ok, note = _validate_expectation(case.get("expect"), name=name, exception=None)
            result["status"] = "pass" if ok else "fail"
            result["reason"] = f"worker_error:{build_reason}; {note}"
            result["name"] = name
        elif bool(case.get("skip_if_no_rdkit", False)):
            result["status"] = "skip"
            result["reason"] = f"worker_error:{build_reason}"
        else:
            result["status"] = "error"
            result["reason"] = f"worker_error:{build_reason}"
    elif build_status is not None:
        result["status"] = "error"
        result["reason"] = str(build_reason or build_status)
    else:
        name_start = time.perf_counter()
        try:
            name = iupac_name(graph, opts)
        except Exception as exc:  # pragma: no cover - defensivo
            exception = exc
            name = None
        name_ms = (time.perf_counter() - name_start) * 1000.0
        ok, note = _validate_expectation(case.get("expect"), name=name, exception=exception)
        result["name"] = str(name or "")
        if ok:
            result["status"] = "pass"
            result["reason"] = note
        else:
            result["status"] = "fail"
            if exception is not None:
                result["reason"] = f"{type(exception).__name__}: {exception}"
            else:
                result["reason"] = note
        if exception is not None:
            result["exception_type"] = type(exception).__name__
            result["exception"] = str(exception)
            if isinstance(exception, ChemNameNotSupported):
                result["not_supported"] = True

    duration_ms = (time.perf_counter() - started) * 1000.0
    result["build_ms"] = round(build_ms, 3)
    result["name_ms"] = round(name_ms, 3)
    result["duration_ms"] = round(duration_ms, 3)

    max_ms = case.get("max_ms")
    if max_ms is not None and result["status"] == "pass":
        try:
            if float(duration_ms) > float(max_ms):
                result["status"] = "fail"
                result["reason"] = f"max_ms_exceeded:{duration_ms:.2f}>{float(max_ms):.2f}"
        except Exception:
            pass
    return result


def _write_junit(results: list[dict[str, Any]], output_path: Path) -> None:
    failures = sum(1 for r in results if r.get("status") == "fail")
    errors = sum(1 for r in results if r.get("status") == "error")
    skipped = sum(1 for r in results if r.get("status") == "skip")
    total_time = sum(float(r.get("duration_ms", 0.0)) for r in results) / 1000.0

    suite = ET.Element(
        "testsuite",
        {
            "name": "chemname_acceptance",
            "tests": str(len(results)),
            "failures": str(failures),
            "errors": str(errors),
            "skipped": str(skipped),
            "time": f"{total_time:.6f}",
        },
    )
    for item in results:
        case = ET.SubElement(
            suite,
            "testcase",
            {
                "classname": "chemname.acceptance",
                "name": str(item.get("id", "")),
                "time": f"{float(item.get('duration_ms', 0.0)) / 1000.0:.6f}",
            },
        )
        status = str(item.get("status", "error"))
        reason = str(item.get("reason", ""))
        if status == "skip":
            ET.SubElement(case, "skipped", {"message": reason})
        elif status == "fail":
            node = ET.SubElement(case, "failure", {"message": reason})
            node.text = json.dumps(item, ensure_ascii=False, indent=2)
        elif status == "error":
            node = ET.SubElement(case, "error", {"message": reason})
            node.text = json.dumps(item, ensure_ascii=False, indent=2)
    tree = ET.ElementTree(suite)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)


def _print_table(results: list[dict[str, Any]], top_slowest: int = 10) -> None:
    print("id                           status   total_ms  name")
    print("-" * 96)
    for item in results:
        name = str(item.get("name", ""))
        if len(name) > 52:
            name = name[:49] + "..."
        print(
            f"{str(item.get('id', ''))[:28]:<28} "
            f"{str(item.get('status', '')):<7} "
            f"{float(item.get('duration_ms', 0.0)):>8.2f}  "
            f"{name}"
        )
    print("")
    ranked = sorted(results, key=lambda row: float(row.get("duration_ms", 0.0)), reverse=True)
    print(f"Top {max(1, int(top_slowest))} más lentos:")
    for item in ranked[: max(1, int(top_slowest))]:
        print(
            f"  - {item.get('id')}: {float(item.get('duration_ms', 0.0)):.2f} ms "
            f"(status={item.get('status')}, reason={item.get('reason')})"
        )


def run_acceptance(
    *,
    cases_path: str | Path = DEFAULT_CASES_PATH,
    report_path: str | Path = DEFAULT_REPORT_PATH,
    junit_path: str | Path | None = None,
    top_slowest: int = 10,
    only_ids: set[str] | None = None,
    limit: int | None = None,
    quiet: bool = False,
) -> dict[str, Any]:
    """Ejecuta todos los casos y genera reportes."""
    cases = load_cases(cases_path)
    if only_ids:
        cases = [case for case in cases if str(case.get("id", "")) in only_ids]
    if limit is not None and limit > 0:
        cases = cases[: int(limit)]

    rdkit_available = is_rdkit_worker_available()
    results: list[dict[str, Any]] = []
    for case in cases:
        results.append(
            evaluate_case(
                case,
                cases_path,
                rdkit_available=rdkit_available,
            )
        )

    summary = {
        "cases_path": str(Path(cases_path).resolve()),
        "report_path": str(Path(report_path).resolve()),
        "rdkit_worker_available": bool(rdkit_available),
        "total": len(results),
        "passed": sum(1 for row in results if row.get("status") == "pass"),
        "failed": sum(1 for row in results if row.get("status") == "fail"),
        "skipped": sum(1 for row in results if row.get("status") == "skip"),
        "errors": sum(1 for row in results if row.get("status") == "error"),
        "slowest": sorted(
            (
                {
                    "id": row.get("id"),
                    "duration_ms": row.get("duration_ms"),
                    "status": row.get("status"),
                }
                for row in results
            ),
            key=lambda row: float(row.get("duration_ms", 0.0)),
            reverse=True,
        )[: max(1, int(top_slowest))],
        "results": results,
    }

    report_file = Path(report_path)
    report_file.parent.mkdir(parents=True, exist_ok=True)
    report_file.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    if junit_path:
        _write_junit(results, Path(junit_path))

    if not quiet:
        _print_table(results, top_slowest=top_slowest)
        print(
            "Summary: "
            f"passed={summary['passed']} failed={summary['failed']} "
            f"skipped={summary['skipped']} errors={summary['errors']} "
            f"(rdkit_worker_available={summary['rdkit_worker_available']})"
        )
        print(f"Report JSON: {report_file}")
        if junit_path:
            print(f"JUnit XML: {Path(junit_path)}")
    return summary


def _parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Acceptance tests para ChemName")
    parser.add_argument(
        "--cases",
        default=str(DEFAULT_CASES_PATH),
        help="Ruta YAML/JSON de casos (default: tests/data/chemname_acceptance_cases.yml)",
    )
    parser.add_argument(
        "--report",
        default=str(DEFAULT_REPORT_PATH),
        help="Ruta de salida JSON (default: artifacts/chemname_acceptance_report.json)",
    )
    parser.add_argument(
        "--junit",
        default=None,
        help="Ruta opcional de salida JUnit XML",
    )
    parser.add_argument("--top-slowest", type=int, default=10, help="Cantidad de casos lentos a mostrar")
    parser.add_argument("--limit", type=int, default=None, help="Limita cantidad de casos ejecutados")
    parser.add_argument(
        "--only",
        default="",
        help="IDs separados por coma para ejecutar subset",
    )
    parser.add_argument("--quiet", action="store_true", help="No imprimir tabla de resultados")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    only_ids = {token.strip() for token in str(args.only).split(",") if token.strip()} or None
    try:
        summary = run_acceptance(
            cases_path=args.cases,
            report_path=args.report,
            junit_path=args.junit,
            top_slowest=max(1, int(args.top_slowest)),
            only_ids=only_ids,
            limit=args.limit,
            quiet=bool(args.quiet),
        )
    except Exception as exc:  # pragma: no cover - CLI defensivo
        print(f"[chemname_acceptance] ERROR: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 2
    if int(summary.get("failed", 0)) or int(summary.get("errors", 0)):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
