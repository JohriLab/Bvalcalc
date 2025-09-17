"""
Header utilities for Bvalcalc output files.

Simple header generation for output files with version, command, warnings, and format.
"""

import sys
from typing import List, Optional
from dataclasses import dataclass


@dataclass
class HeaderInfo:
    """Container for header information."""
    file_type: str
    command: Optional[str] = None
    warnings: List[str] = None
    data_format: Optional[str] = None
    
    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


def generate_header(info: HeaderInfo) -> List[str]:
    """
    Generate header lines from HeaderInfo object.
    
    Args:
        info: HeaderInfo object with information to include
        
    Returns:
        List of header lines (with # prefix)
    """
    lines = []
    
    # File identification
    try:
        from Bvalcalc import __version__
        lines.append(f"# Bvalcalc v{__version__}")
    except ImportError:
        lines.append("# Bvalcalc v(version not found)")
    
    # Command information
    if info.command:
        # Remove the script path from the beginning since it's already in the first line
        command = info.command
        # Check for various ways the script might be called
        if "Bvalcalc/__main__.py" in command:
            # Extract just the arguments after the script path
            parts = command.split()
            if len(parts) > 1:
                command = " ".join(parts[1:])
        elif command.startswith("python -m Bvalcalc "):
            command = command[19:]  # Remove "python -m Bvalcalc "
        elif command.startswith("bvalcalc "):
            command = command[9:]  # Remove "bvalcalc "
        lines.append(f"# {command}")
    
    # Warnings
    if info.warnings:
        for warning in info.warnings:
            lines.append(f"# WARNING: {warning}")
    
    # Data format
    if info.data_format:
        lines.append(f"# Format: {info.data_format}")
    
    return lines


def write_headers_to_file(file_path: str, header_lines: List[str], mode: str = 'w'):
    """
    Write header lines to the specified file.
    
    Args:
        file_path: Path to the file
        header_lines: List of header lines (with # prefix)
        mode: File open mode ('w' or 'a')
    """
    with open(file_path, mode) as f:
        for line in header_lines:
            f.write(line + "\n")


def create_header_info_from_args(args, file_type: str, description: str = None) -> HeaderInfo:
    """
    Create HeaderInfo from command line arguments.
    
    Args:
        args: Command line arguments object
        file_type: Type of file being generated
        description: Optional description
        
    Returns:
        HeaderInfo object with information from args
    """
    info = HeaderInfo(
        file_type=file_type,
        data_format="Chromosome,Start,B"
    )
    
    # Build command string
    command_parts = [sys.argv[0]]
    for arg in sys.argv[1:]:
        if ' ' in arg:
            command_parts.append(f'"{arg}"')
        else:
            command_parts.append(arg)
    info.command = ' '.join(command_parts)
    
    return info