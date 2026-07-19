import yaml
import pytest

def test_exceptions_no_growth():
    with open('architecture/modules.yml', 'r') as f:
        data = yaml.safe_load(f)
    
    modules = data['modules']
    total_exceptions = 0
    for m in modules:
        total_exceptions += len(m.get('temporary_exceptions', []))
    
    # For the purpose of this task, we expect a certain number of exceptions
    # based on our current implementation.
    # Let's count them manually for now to establish the baseline.
    # M01: 6
    # M03: 2
    # M04: 7
    # M15: 2
    # Total = 17
    
    expected_count = 17
    assert total_exceptions == expected_count, f"Unexpected number of exceptions: found {total_exceptions}, expected {expected_count}"

if __name__ == "__main__":
    test_exceptions_no_growth()
