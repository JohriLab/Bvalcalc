"""
Bvalcalc: calculate relative diversity (B) under background selection.
"""

__version__ = "0.2.0"

# Expose main entry point
from .cli import main

__all__ = [
    "main",
    "__version__",
]