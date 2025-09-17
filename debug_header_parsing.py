#!/usr/bin/env python3
"""
Debug script for header parsing.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'Bvalcalc'))

from Bvalcalc.utils.header_utils import parse_headers, extract_header_info

# Test content
test_content = """# Bvalcalc Header Format v1.0
# Generated: 2024-01-15 14:30:25
# File Type: BED
# Description: Test BED file with headers
# Command: bvalcalc --genome --params test_params.py --bedgff test.bed
#
# WARNING: This is a test warning
#
# Input Files:
#   - test.bed: Test BED file
#   - test_params.py: Test parameters
#
# Parameters:
#   - r: 1.2e-8
#   - g: 1.0e-6
#
# Data: Chromosome,Start,End
chr2L,1000,2000
chr2L,3000,4000
"""

def parse_headers_from_content(content):
    """Parse headers from content string."""
    lines = content.strip().split('\n')
    header_lines = []
    data_start = 0
    
    for line_num, line in enumerate(lines, 1):
        line = line.strip()
        if line.startswith('#'):
            header_lines.append(line)
        else:
            data_start = line_num
            break
    
    return header_lines, data_start

print("Header lines:")
header_lines, data_start = parse_headers_from_content(test_content)
for i, line in enumerate(header_lines):
    print(f"{i:2d}: {line}")
    if "r:" in line:
        print(f"    -> Parameter line: '{line}' (starts with '  - '? {line.startswith('#  - ')})")

print(f"\nData starts at line: {data_start}")

print("\nParsing header info...")
# Add debug output to the parsing function
def debug_extract_header_info(header_lines):
    from Bvalcalc.utils.header_utils import HeaderInfo
    info = HeaderInfo(file_type="unknown")
    
    in_parameters_section = False
    in_input_files_section = False
    
    for line in header_lines:
        original_line = line
        line = line[1:].strip()  # Remove # and strip whitespace
        
        if "r:" in line or "g:" in line:
            print(f"  -> Processing line: '{original_line}' -> '{line}'")
            print(f"     starts with '  - '? {line.startswith('  - ')}")
            print(f"     starts with '   - '? {line.startswith('   - ')}")
            print(f"     in_parameters_section? {in_parameters_section}")
            print(f"     contains ':'? {':' in line}")
        
        if line.startswith("Parameters:"):
            in_parameters_section = True
            in_input_files_section = False
            print(f"  -> Entering parameters section")
            continue
        elif line.startswith("Input Files:"):
            in_input_files_section = True
            in_parameters_section = False
            print(f"  -> Entering input files section")
            continue
        elif line.startswith("- ") and in_parameters_section and ":" in line:
            print(f"  -> Found parameter line: '{line}' (in_parameters_section: {in_parameters_section})")
            # Parameter entry (indented with spaces)
            parts = line[2:].split(":", 1)  # Remove "- " prefix
            if len(parts) == 2:
                param_name = parts[0].strip()
                param_value = parts[1].strip()
                info.parameters[param_name] = param_value
                print(f"  -> Added parameter: {param_name} = {param_value}")
        # Don't reset section flags on empty lines - they're just formatting
    
    return info

header_info = debug_extract_header_info(header_lines)

print(f"File type: {header_info.file_type}")
print(f"Description: {header_info.description}")
print(f"Warnings: {header_info.warnings}")
print(f"Input files: {header_info.input_files}")
print(f"Parameters: {header_info.parameters}")
