"""
bvalcalc: calculate relative diversity (B) under background selection.
"""

__version__ = "1.3.0"

# Expose main entry point
from .cli import main
from .core.calculateB import calculateB_linear, calculateB_recmap, calculateB_unlinked, calculateB_hri, get_params


__all__ = [
    "get_params", "calculateB_linear", "calculateB_unlinked", "calculateB_hri",
    "main",
    "__version__",
]