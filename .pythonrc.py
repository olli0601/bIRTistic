# Python startup script for bIRTistic
# Automatically adds the python module to sys.path

import sys
from pathlib import Path
import pandas as pd

# Auto-detect project root
cwd = Path.cwd()
if cwd.name == 'scripts-py':
    project_root = cwd.parent
elif (cwd / 'scripts-py').exists():
    project_root = cwd
elif (cwd / 'python').exists():
    project_root = cwd
else:
    # Try to find git root
    current = cwd
    while current != current.parent:
        if (current / '.git').exists():
            project_root = current
            break
        current = current.parent
    else:
        project_root = cwd

# Add python module to path
python_path = str(project_root / 'python')
if python_path not in sys.path:
    sys.path.insert(0, python_path)

# Global pandas display options for interactive sessions
pd.set_option('display.max_columns', None)

print(f"✓ bIRTistic Python environment ready")
print(f"  Project root: {project_root}")
print(f"  Python path: {python_path}")
print(f"  pandas display.max_columns: {pd.get_option('display.max_columns')}")
print(f"  Available modules: data_loading, model_fitting, analysis")
print()
