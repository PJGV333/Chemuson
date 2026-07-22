from collections import Counter
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
import re
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CATALOG_PATH = REPO_ROOT / "architecture" / "modules.yml"
@dataclass(frozen=True, order=True)
class ExceptionIdentity:
    source_id: str; target_id: str; file: str; import_path: str; type_checking_only: bool
def _identity(source_id, target_id, file, import_path, type_checking_only):
    return ExceptionIdentity(source_id,target_id,PurePosixPath(file.strip().replace("\\","/")).as_posix(),re.sub(r"\s+"," ",import_path).strip(),type_checking_only)
PERSISTENCE_EXCEPTION = _identity("M01","M09","src/chemuson/chemio/persistence.py","from chemuson.gui.canvas import ChemusonCanvas",True)
FROZEN_EXCEPTION_BASELINE_ROWS = (PERSISTENCE_EXCEPTION,)
FROZEN_EXCEPTION_BASELINE = frozenset(FROZEN_EXCEPTION_BASELINE_ROWS)
ELIMINATED_M01_M02_EXCEPTIONS = (
 _identity("M01","M02","src/chemuson/chemio/depiction_candidates.py","from chemuson.clean2d.geometry import count_crossings, cycle_basis, segments_intersect",False),
 _identity("M01","M02","src/chemuson/chemio/depiction_candidates.py","from chemuson.clean2d.safety import has_cycles, min_nonbonded_distance, ring_degeneracy_score",False),
 _identity("M01","M02","src/chemuson/chemio/rdkit_io.py","from chemuson.clean2d.geometry import apply_coords_in_place",False),
 _identity("M01","M02","src/chemuson/chemio/rdkit_io.py","from chemuson.clean2d.scaffold_depiction import scaffold_depiction_candidates",False),
 _identity("M01","M02","src/chemuson/chemio/rdkit_io.py","from chemuson.clean2d.block_unwrap import block_unwrap_layout",False),
)
def audit(rows):
    current=tuple(rows); unexpected=tuple(sorted(set(current)-FROZEN_EXCEPTION_BASELINE)); growth=Counter(x.source_id for x in current); return unexpected, growth
def test_frozen_baseline_is_exact():
 assert FROZEN_EXCEPTION_BASELINE_ROWS == (PERSISTENCE_EXCEPTION,)
def test_real_catalog_matches_frozen_exception_baseline():
 rows=tuple(_identity(e["source_id"],e["target_id"],e["file"],e["import_path"],e["type_checking_only"]) for m in yaml.safe_load(CATALOG_PATH.read_text())["modules"] for e in m["temporary_exceptions"] )
 assert set(rows) <= FROZEN_EXCEPTION_BASELINE
 assert len(rows)==1 and rows[0] == PERSISTENCE_EXCEPTION
def test_eliminated_m01_m02_exceptions_cannot_reappear():
 unexpected,growth=audit(ELIMINATED_M01_M02_EXCEPTIONS)
 assert unexpected == tuple(sorted(ELIMINATED_M01_M02_EXCEPTIONS)) and growth["M01"] == 5
def test_partial_and_replacement_m01_m02_are_rejected():
 unexpected,_=audit((ELIMINATED_M01_M02_EXCEPTIONS[0],)); assert unexpected
 unexpected,growth=audit((ELIMINATED_M01_M02_EXCEPTIONS[0],)); assert unexpected and growth["M01"] == 1
def test_removing_persistence_exception_is_allowed():
 unexpected,growth=audit(()); assert not unexpected and not growth
def test_normalized_historical_identity_is_still_rejected():
 row=ELIMINATED_M01_M02_EXCEPTIONS[0]; unexpected,_=audit((_identity(row.source_id,row.target_id,r"src\chemuson\chemio\depiction_candidates.py","from chemuson.clean2d.geometry\n import count_crossings, cycle_basis, segments_intersect",False),)); assert unexpected == (row,)
