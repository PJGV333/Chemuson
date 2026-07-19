import os
import pytest

def test_no_tools_in_src():
    for root, _, files in os.walk('src/chemuson/'):
        for file in files:
            if file.endswith('.py'):
                full_path = os.path.join(root, file)
                with open(full_path, 'r') as f:
                    content = f.read()
                    if 'import tools' in content or 'from tools' in content:
                        pytest.fail(f"Architecture violation: {full_path} imports 'tools'")
                    if 'import chemuson.tools' in content or 'from chemuson.tools' in content:
                        pytest.fail(f"Architecture violation: {full_path} imports 'chemuson.tools'")

if __name__ == "__main__":
    test_no_tools_in_src()
