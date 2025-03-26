import os
import sys

def pytest_sessionstart(session):
    # Ensure project root is in sys.path
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    # Set env var for pop params
    if "BCALC_POP_PARAMS" not in os.environ:
        param_path = os.path.join(root_dir, "ExampleParams.py")
        if not os.path.exists(param_path):
            raise FileNotFoundError(f"Could not find ExampleParams.py at expected path: {param_path}")
        os.environ["BCALC_POP_PARAMS"] = param_path
