#!/usr/bin/env python3
"""
Test script for header functionality in Bvalcalc.

This script tests the new header parsing and generation functionality
across different file types.
"""

import os
import sys
import tempfile
import numpy as np

# Add the Bvalcalc package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'Bvalcalc'))

from Bvalcalc.utils.header_utils import (
    HeaderInfo, generate_header, parse_headers, extract_header_info,
    create_header_info_from_args, add_warning_to_headers
)
from Bvalcalc.utils.load_bed_gff import get_bed_gff_headers, get_bed_gff_header_info
from Bvalcalc.utils.load_rec_map import get_rec_map_headers, get_rec_map_header_info
from Bvalcalc.utils.load_chr_sizes import get_chr_sizes_headers, get_chr_sizes_header_info


def test_header_generation():
    """Test header generation functionality."""
    print("Testing header generation...")
    
    # Create a test HeaderInfo object
    info = HeaderInfo(
        file_type="B-map",
        description="Test B-value calculations",
        mode="genome",
        command="bvalcalc --genome --params test_params.py --bedgff test.bed --out test.csv",
        warnings=["Test warning message"],
        input_files=[("test.bed", "Test BED file"), ("test_params.py", "Test parameters")],
        parameters={"r": "1.2e-8", "g": "1.0e-6", "Nanc": "1000000"},
        data_format="Chromosome,Start,B"
    )
    
    # Generate header lines
    header_lines = generate_header(info)
    
    print("Generated header:")
    for line in header_lines:
        print(f"  {line}")
    
    # Verify header format
    assert header_lines[0] == "# Bvalcalc v1.0.0"
    assert any("WARNING: Test warning message" in line for line in header_lines)
    assert any("Format: Chromosome,Start,B" in line for line in header_lines)
    
    print("✓ Header generation test passed!")


def test_header_parsing():
    """Test header parsing functionality."""
    print("\nTesting header parsing...")
    
    # Create a test file with headers
    test_content = """# Bvalcalc v1.0.0
# --genome --params test_params.py --bedgff test.bed
# WARNING: This is a test warning
# Format: Chromosome,Start,End
chr2L,1000,2000
chr2L,3000,4000
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(test_content)
        temp_file = f.name
    
    try:
        # Test header parsing
        header_lines, data_start = parse_headers(temp_file)
        
        print(f"Found {len(header_lines)} header lines")
        print(f"Data starts at line {data_start}")
        
        # Test header info extraction
        header_info = extract_header_info(header_lines)
        
        print(f"File type: {header_info.file_type}")
        print(f"Description: {header_info.description}")
        print(f"Warnings: {header_info.warnings}")
        print(f"Input files: {header_info.input_files}")
        print(f"Parameters: {header_info.parameters}")
        
        # Verify parsed information
        assert header_info.file_type == "bvalcalc"
        assert "This is a test warning" in header_info.warnings
        assert header_info.data_format == "Chromosome,Start,End"
        
        print("✓ Header parsing test passed!")
        
    finally:
        os.unlink(temp_file)


def test_file_parsing_with_headers():
    """Test that file parsing functions handle headers correctly."""
    print("\nTesting file parsing with headers...")
    
    # Test BED file with headers
    bed_content = """# Bvalcalc v1.0.0
# Format: Chromosome,Start,End
chr2L,1000,2000
chr2L,3000,4000
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        f.write(bed_content)
        bed_file = f.name
    
    try:
        # Test BED header functions
        headers = get_bed_gff_headers(bed_file)
        header_info = get_bed_gff_header_info(bed_file)
        
        print(f"BED headers: {len(headers)} lines")
        print(f"BED file type: {header_info.file_type}")
        
        assert len(headers) > 0
        assert header_info.file_type == "bvalcalc"
        
        print("✓ BED file parsing with headers test passed!")
        
    finally:
        os.unlink(bed_file)
    
    # Test CSV file with headers
    csv_content = """# Bvalcalc v1.0.0
# Format: Chromosome,Start,Rate
chr2L,1000,1.2e-8
chr2L,2000,1.5e-8
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write(csv_content)
        csv_file = f.name
    
    try:
        # Test CSV header functions
        headers = get_rec_map_headers(csv_file)
        header_info = get_rec_map_header_info(csv_file)
        
        print(f"CSV headers: {len(headers)} lines")
        print(f"CSV file type: {header_info.file_type}")
        
        assert len(headers) > 0
        assert header_info.file_type == "bvalcalc"
        
        print("✓ CSV file parsing with headers test passed!")
        
    finally:
        os.unlink(csv_file)


def test_warning_handling():
    """Test warning handling in headers."""
    print("\nTesting warning handling...")
    
    # Create header with warnings
    header_lines = [
        "# Bvalcalc v1.0.0",
        "# WARNING: First warning",
        "# WARNING: Second warning"
    ]
    
    # Add a new warning
    new_header_lines = add_warning_to_headers(header_lines, "New warning message")
    
    print("Updated header with new warning:")
    for line in new_header_lines:
        print(f"  {line}")
    
    # Verify warning was added
    warning_lines = [line for line in new_header_lines if "WARNING:" in line]
    assert len(warning_lines) == 3  # Original 2 + new 1
    assert any("New warning message" in line for line in warning_lines)
    
    print("✓ Warning handling test passed!")


def main():
    """Run all tests."""
    print("Running header functionality tests...\n")
    
    try:
        test_header_generation()
        test_header_parsing()
        test_file_parsing_with_headers()
        test_warning_handling()
        
        print("\n🎉 All tests passed! Header functionality is working correctly.")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
