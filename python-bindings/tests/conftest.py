from pathlib import Path
import sys


PYTHON_BINDINGS_DIRECTORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_BINDINGS_DIRECTORY / "build"))
