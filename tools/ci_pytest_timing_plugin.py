"""Write completed pytest phase durations as newline-delimited JSON."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any


class ShardTimingPlugin:
    """Persist completed node timings incrementally so timeouts retain evidence."""

    def __init__(self, output_path: Path) -> None:
        self.output_path = output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("", encoding="utf-8")
        self._pending: dict[str, dict[str, Any]] = {}
        self._emitted: set[str] = set()

    def pytest_runtest_logreport(self, report: Any) -> None:
        nodeid = str(report.nodeid)
        if nodeid in self._emitted:
            return
        entry = self._pending.setdefault(
            nodeid,
            {"nodeid": nodeid, "test_file": nodeid.split("::", 1)[0], "phases": {}, "outcome": "passed"},
        )
        entry["phases"][str(report.when)] = max(0.0, float(report.duration))
        if report.failed:
            entry["outcome"] = "failed"
        elif report.skipped and entry["outcome"] == "passed":
            entry["outcome"] = "skipped"

        if report.when == "teardown" or (
            report.when == "setup" and (report.failed or report.skipped)
        ):
            entry["duration_seconds"] = round(sum(entry["phases"].values()), 6)
            with self.output_path.open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(entry, sort_keys=True) + "\n")
                stream.flush()
            del self._pending[nodeid]
            self._emitted.add(nodeid)


def pytest_configure(config: Any) -> None:
    output = os.environ.get("CHEMUSON_PYTEST_TIMING_PATH")
    if output:
        config.pluginmanager.register(
            ShardTimingPlugin(Path(output)), "chemuson-shard-timing"
        )
