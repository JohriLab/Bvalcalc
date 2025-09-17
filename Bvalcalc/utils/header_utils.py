"""
Header utilities for Bvalcalc file I/O.

This module provides functionality to parse, generate, and manage headers
in Bvalcalc input and output files. Headers use # comments to provide
metadata, warnings, and input commands.
"""

import os
import datetime
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass


@dataclass
class HeaderInfo:
    """Container for header information."""
    file_type: str
    generated: Optional[str] = None
    command: Optional[str] = None
    mode: Optional[str] = None
    warnings: List[str] = None
    input_files: List[Tuple[str, str]] = None  # (file_path, description)
    parameters: Dict[str, str] = None
    data_format: Optional[str] = None
    description: Optional[str] = None
    
    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []
        if self.input_files is None:
            self.input_files = []
        if self.parameters is None:
            self.parameters = {}


def parse_headers(file_path: str) -> Tuple[List[str], int]:
    """
    Parse headers from a file and return header lines and data start line.
    
    Args:
        file_path: Path to the file to parse
        
    Returns:
        Tuple of (header_lines, data_start_line_number)
    """
    header_lines = []
    data_start_line = 0
    
    with open(file_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line.startswith('#'):
                header_lines.append(line)
            else:
                data_start_line = line_num
                break
    
    return header_lines, data_start_line


def extract_header_info(header_lines: List[str]) -> HeaderInfo:
    """
    Extract structured information from header lines.
    
    Args:
        header_lines: List of header lines (starting with #)
        
    Returns:
        HeaderInfo object with extracted information
    """
    info = HeaderInfo(file_type="unknown")
    
    in_parameters_section = False
    in_input_files_section = False
    
    for line in header_lines:
        line = line[1:].strip()  # Remove # and strip whitespace
        
        if line.startswith("Bvalcalc Header Format"):
            info.file_type = "bvalcalc"
        elif line.startswith("Generated:"):
            info.generated = line.split(":", 1)[1].strip()
        elif line.startswith("Command:"):
            info.command = line.split(":", 1)[1].strip()
        elif line.startswith("Mode:"):
            info.mode = line.split(":", 1)[1].strip()
        elif line.startswith("WARNING:"):
            warning = line.split(":", 1)[1].strip()
            info.warnings.append(warning)
        elif line.startswith("File Type:"):
            info.file_type = line.split(":", 1)[1].strip()
        elif line.startswith("Description:"):
            info.description = line.split(":", 1)[1].strip()
        elif line.startswith("Data:"):
            info.data_format = line.split(":", 1)[1].strip()
        elif line.startswith("Input Files:"):
            in_input_files_section = True
            in_parameters_section = False
            continue
        elif line.startswith("Parameters:"):
            in_parameters_section = True
            in_input_files_section = False
            continue
        elif line.startswith("- ") and in_input_files_section:
            # Input file entry
            parts = line[2:].split(":", 1)
            if len(parts) == 2:
                file_path = parts[0].strip()
                description = parts[1].strip()
                info.input_files.append((file_path, description))
        elif line.startswith("- ") and in_parameters_section and ":" in line:
            # Parameter entry (indented with spaces)
            parts = line[2:].split(":", 1)  # Remove "- " prefix
            if len(parts) == 2:
                param_name = parts[0].strip()
                param_value = parts[1].strip()
                info.parameters[param_name] = param_value
        # Don't reset section flags on empty lines - they're just formatting
    
    return info


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
    lines.append("# Bvalcalc Header Format v1.0")
    if info.generated:
        lines.append(f"# Generated: {info.generated}")
    else:
        lines.append(f"# Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if info.file_type:
        lines.append(f"# File Type: {info.file_type}")
    
    if info.description:
        lines.append(f"# Description: {info.description}")
    
    # Command information
    if info.command:
        lines.append(f"# Command: {info.command}")
    
    if info.mode:
        lines.append(f"# Mode: {info.mode}")
    
    # Warnings
    if info.warnings:
        lines.append("")
        for warning in info.warnings:
            lines.append(f"# WARNING: {warning}")
    
    # Input files
    if info.input_files:
        lines.append("")
        lines.append("# Input Files:")
        for file_path, description in info.input_files:
            lines.append(f"#   - {file_path}: {description}")
    
    # Parameters
    if info.parameters:
        lines.append("")
        lines.append("# Parameters:")
        for param_name, param_value in info.parameters.items():
            lines.append(f"#   - {param_name}: {param_value}")
    
    # Data format
    if info.data_format:
        lines.append("")
        lines.append(f"# Data: {info.data_format}")
    
    # Add empty line before data
    lines.append("")
    
    return lines


def write_headers_to_file(file_path: str, header_lines: List[str], mode: str = 'w'):
    """
    Write header lines to a file.
    
    Args:
        file_path: Path to the file to write
        header_lines: List of header lines (with # prefix)
        mode: File open mode ('w' for write, 'a' for append)
    """
    with open(file_path, mode) as f:
        for line in header_lines:
            f.write(line + '\n')


def preserve_headers_from_input(input_file: str) -> List[str]:
    """
    Extract and preserve headers from an input file.
    
    Args:
        input_file: Path to the input file
        
    Returns:
        List of header lines (with # prefix)
    """
    if not os.path.exists(input_file):
        return []
    
    header_lines, _ = parse_headers(input_file)
    return header_lines


def add_warning_to_headers(header_lines: List[str], warning: str) -> List[str]:
    """
    Add a warning to existing header lines.
    
    Args:
        header_lines: Existing header lines
        warning: Warning message to add
        
    Returns:
        Updated header lines with warning added
    """
    # Find the best place to insert the warning
    warning_line = f"# WARNING: {warning}"
    
    # If there are already warnings, add after the last one
    last_warning_idx = -1
    for i, line in enumerate(header_lines):
        if line.startswith("# WARNING:"):
            last_warning_idx = i
    
    if last_warning_idx >= 0:
        # Insert after the last warning
        header_lines.insert(last_warning_idx + 1, warning_line)
    else:
        # Find a good place to insert warnings section
        # Look for empty line or end of file identification section
        insert_idx = 0
        for i, line in enumerate(header_lines):
            if line.startswith("# Generated:") or line.startswith("# File Type:"):
                insert_idx = i + 1
            elif line.strip() == "#" and i > insert_idx:
                insert_idx = i
                break
        
        if insert_idx < len(header_lines):
            header_lines.insert(insert_idx, "")
            header_lines.insert(insert_idx + 1, warning_line)
        else:
            header_lines.append("")
            header_lines.append(warning_line)
    
    return header_lines


def create_header_info_from_args(args, file_type: str, description: str = None) -> HeaderInfo:
    """
    Create HeaderInfo from command line arguments.
    
    Args:
        args: Parsed command line arguments
        file_type: Type of file being created
        description: Optional description of the file
        
    Returns:
        HeaderInfo object with information from args
    """
    info = HeaderInfo(
        file_type=file_type,
        description=description,
        mode=getattr(args, 'mode', None)
    )
    
    # Build command string
    import sys
    command_parts = [sys.argv[0]]
    for arg in sys.argv[1:]:
        if ' ' in arg:
            command_parts.append(f'"{arg}"')
        else:
            command_parts.append(arg)
    info.command = ' '.join(command_parts)
    
    # Add input files
    if hasattr(args, 'bedgff') and args.bedgff:
        info.input_files.append((args.bedgff, "BED/GFF annotations"))
    
    if hasattr(args, 'rec_map') and args.rec_map:
        info.input_files.append((args.rec_map, "Recombination map"))
    
    if hasattr(args, 'gc_map') and args.gc_map:
        info.input_files.append((args.gc_map, "Gene conversion map"))
    
    if hasattr(args, 'chr_sizes') and args.chr_sizes:
        info.input_files.append((args.chr_sizes, "Chromosome sizes"))
    
    if hasattr(args, 'params') and args.params:
        info.input_files.append((args.params, "Parameters file"))
    
    return info
