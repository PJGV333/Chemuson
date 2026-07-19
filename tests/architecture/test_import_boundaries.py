import ast
import os
import yaml
import pytest
from pathlib import Path

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
                # Path is like src/chemuson/gui/
                rel_path = p.relative_to('src/chemuson/')
                if str(rel_path) == '.':
                    package_prefix = 'chemuson'
                else:
                    package_prefix = f"chemuson.{rel_path.as_posix().replace('/', '.')}".rstrip('.')
                mapping[package_prefix] = m_id
            elif p.is_file():
                # Path is like src/chemuson/gui/__init__.py or src/chemuson/core/model.py
                p_rel = p.relative_to('src')
                parts = list(p_rel.parts)
                if parts[-1] == '__init__.py':
                    parts.pop() # remove __init__.py from parts
                else:
                    parts[-1] = p_rel.stem # remove .py
                package_prefix = '.'.join(parts)
                mapping[package_prefix] = m_id
    return mapping

def get_exceptions_map():
    modules = load_catalog()
    exceptions = []
    for m in modules:
        for exc in m.get('temporary_exceptions', []):
            exceptions.append(exc)
    return exceptions

def get_forbidden_for_module(m_id):
    modules = load_catalog()
    for m in modules:
        if m['id'] == m_id:
            return set(m.get('forbidden_dependencies', []))
    return set()

def get_target_dependencies(m_id):
    modules = load_catalog()
    for m in modules:
        if m['id'] == m_id:
            return set(m.get('target_dependencies', []))
    return set()

def get_module_id_from_pkg(pkg_name, mapping):
    parts = pkg_name.split('.')
    for i in range(len(parts), 0, -1):
        prefix = '.'.join(parts[:i])
        if prefix in mapping:
            return mapping[prefix]
    return None

class ImprovedImportVisitor(ast.NodeVisitor):
    def __init__(self, filename, lines):
        self.filename = filename
        self.lines = lines
        self.imports = [] 

    def visit_ImportFrom(self, node):
        if node.module and node.module.startswith('chemuson'):
            line = self.lines[node.lineno - 1].strip()
            self.imports.append(('ImportFrom', node.module, line, node.lineno, node))
        self.generic_visit(node)

    def visit_Import(self, node):
        for alias in node.names:
            if alias.name.startswith('chemuson'):
                line = self.lines[node.lineno - 1].strip()
                self.imports.append(('Import', alias.name, line, node.lineno, node))
        self.generic_visit(node)

def get_parents_map(tree):
    parents = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parents[child] = node
    return parents

def test_import_boundaries_v3():
    mapping = get_module_map()
    exceptions = get_exceptions_map()
    
    for root, _, files in os.walk('src/chemuson/'):
        for file in files:
            if file.endswith('.py'):
                full_path = os.path.join(root, file)
                rel_path = os.path.relpath(full_path, 'src')
                p_rel = Path(rel_path)
                parts = list(p_rel.parts)
                parts[-1] = p_rel.stem
                pkg_of_file = '.'.join(parts)
                m_id = get_module_id_from_pkg(pkg_of_file, mapping)
                if not m_id:
                    continue
    
                with open(full_path, 'r') as f:
                    content = f.read()
                    lines = content.splitlines()
    
                try:
                    tree = ast.parse(content)
                except Exception as e:
                    pytest.fail(f"Failed to parse {full_path}: {e}")
    
                parents = get_parents_map(tree)
                visitor = ImprovedImportVisitor(full_path, lines)
                visitor.visit(tree)
    
                forbidden = get_forbidden_for_module(m_id)
                target = set(get_target_dependencies(m_id))
    
                for imp_type, imp_pkg, imp_line, line_no, node in visitor.imports:
                    imp_m_id = get_module_id_from_pkg(imp_pkg, mapping)
                    if not imp_m_id:
                        continue
    
                        # Check if it's an exception
                        is_exc = False
                        for exc in exceptions:
                            if exc['source_id'] == m_id and exc['target_id'] == imp_m_id:
                                if exc['import_path'].strip() in imp_line:
                                    if exc.get('type_checking_only', False):
                                        is_tc = False
                                        curr = node
                                        while curr in parents:
                                            curr = parents[curr]
                                            if isinstance(curr, ast.If):
                                                if (isinstance(curr.test, ast.Name) and curr.test.id == 'TYPE_CHECKING') or \
                                                   (isinstance(curr.test, ast.Attribute) and curr.test.attr == 'TYPE_CHECKING'):
                                                    is_tc = True
                                                    break
                                        is_exc = is_tc
                                    else:
                                        is_exc = True
                                    
                                    if is_exc:
                                        break
    
                        if is_exc:
                            continue
    
                        if imp_m_id in forbidden:
                            pytest.fail(f"Architecture violation in {full_path}:{line_no}: Import of forbidden module {imp_m_id}. Statement: `{imp_line}`")
    
                        if imp_m_id != m_id and imp_m_id not in target:
                            print(f"DEBUG: file={full_path}, m_id={m_id}, imp_m_id={imp_m_id}, target={target}")
                            pytest.fail(f"Architecture violation in {full_path}:{line_no}: Import of {imp_m_id} not in target_dependencies. Statement: `{imp_line}`")

test_import_boundaries = test_import_boundaries_v3
