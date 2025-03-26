import os
import sys

def pytest_sessionstart(session):
    # Ensure project root is in sys.path
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if root_dir not in sys.path:
        sys.path.insert(0, root_dir)

    # Set default env var for pop params if not already set
    if "BCALC_POP_PARAMS" not in os.environ:
        default_params = os.path.join(root_dir, "tests", "testparams", "nogcBasicParams.py")
        if not os.path.exists(default_params):
            raise FileNotFoundError(f"Could not find test params at expected path: {default_params}")
        os.environ["BCALC_POP_PARAMS"] = default_params
