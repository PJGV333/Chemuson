import ast
import yaml
from pathlib import Path
import pytest

def load_catalog():
    with open('architecture/modules.yml', 'r') as f:
        return yaml.safe_load(f)['modules']

def get_module_map():
    modules = load_catalog()
    mapping = {}
    for m in modules:
        m_id = m['id']
        for path in m.get('paths', []):
            p = Path(path)
            if p.is_dir():
                rel_path = p.relative_to('src/chemuson/')
                if str(rel_path) == '.':
                    package_prefix = 'chemuson'
                else:
                    package_prefix = f"chemuson.{rel_path.as_posix().replace('/', '.')}".rstrip('.')
                mapping[package_prefix] = m_id
            elif p.is_file():
                p_rel = p.relative_to('src')
                parts = list(p_rel.parts)
                parts[-1] = p_rel.stem
                pkg = '.'.join(parts)
                mapping[pkg] = m_id
    return mapping

def test_public_api_exists():
    modules = load_catalog()
    mapping = get_module_map()
    
    for m in modules:
        public_api = m.get('public_api', [])
        if not public_api:
            continue
            
        # Find the __init__.py for this module
        init_path = None
        for path in m.get('paths', []):
            p = Path(path)
            if p.is_dir():
                init_file = p / "__init__.py"
                if init_file.exists():
                    init_path = init_file
                    break
            elif p.is_file() and p.name == "__init__.py":
                init_path = p
                break
        
        if not init_path:
            # If no __init__.py, maybe it's a single file module?
            # The documentation says "symbols reexported in __init__.py"
            # If it's a single file, maybe the symbols are in the file itself.
            # But the catalog says they are in __init__.py.
            continue

        with open(init_path, 'r') as f:
            content = f.read()
            tree = ast.parse(content)
        
        # Get all defined names in __init__.py
        defined_names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    defined_names.add(alias.name if alias.asname is None else alias.asname)
            elif isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    defined_names.add(alias.name if alias.asname is None else alias.asname)
            elif isinstance(node, ast.FunctionDef):
                defined_names.add(node.name)
            elif isinstance(node, ast.ClassDef):
                defined_names.add(node.name)
        
        for symbol in public_api:
            assert symbol in defined_names, f"Module {m['id']} ({m['name']}) declares '{symbol}' in public_api, but it's not in {init_path}"

