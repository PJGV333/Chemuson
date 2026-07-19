import yaml
import os
import pytest
from pathlib import Path

def load_modules():
    with open('architecture/modules.yml', 'r') as f:
        return yaml.safe_load(f)['modules']

def test_module_ids_range():
    modules = load_modules()
    ids = sorted([m['id'] for m in modules])
    expected = [f"M{i:02d}" for i in range(20)]
    assert ids == expected, f"Expected {expected}, got {ids}"

def test_module_ids_unique():
    modules = load_modules()
    ids = [m['id'] for m in modules]
    assert len(ids) == len(set(ids)), "Module IDs are not unique"

def test_module_paths_exist():
    modules = load_modules()
    for m in modules:
        for path in m.get('paths', []):
            assert Path(path).exists(), f"Path {path} in module {m['id']} does not exist"

def test_module_required_fields_not_empty():
    modules = load_modules()
    required_fields = {'id', 'name', 'title', 'responsibility', 'paths', 'status', 'risk_level'}
    for m in modules:
        for field in required_fields:
            assert field in m, f"Module {m['id']} is missing required field: {field}"
            if field != 'id':
                assert m[field], f"Module {m['id']} has empty field: {field}"

def test_module_status_values():
    modules = load_modules()
    valid_statuses = {'stable', 'evolving', 'legacy', 'empty'}
    for m in modules:
        assert m['status'] in valid_statuses, f"Module {m['id']} has invalid status: {m['status']}"

def test_module_risk_level_values():
    modules = load_modules()
    valid_risks = {'low', 'medium', 'high'}
    for m in modules:
        assert m['risk_level'] in valid_risks, f"Module {m['id']} has invalid risk_level: {m['risk_level']}"
