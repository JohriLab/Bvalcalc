"""
Test header output functionality.
"""

import tempfile
import subprocess
import os
from pathlib import Path


def test_header_output_format():
    """Test that output files have correct header format."""
    
    # Create temporary files for the test
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as bed_file:
        bed_file.write("chr1\t1000\t2000\n")
        bed_file.write("chr1\t3000\t4000\n")
        bed_path = bed_file.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as rec_file:
        # Create a rec map that doesn't start at position 1 (to trigger warning)
        rec_file.write("chr1,5000,1.2e-8\n")
        rec_file.write("chr1,10000,1.5e-8\n")
        rec_path = rec_file.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as chr_file:
        chr_file.write("chr1,20000\n")
        chr_path = chr_file.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as params_file:
        params_file.write("""
# Test parameters
x = 100
Nanc = 1e6/x
Ncur = Nanc
r = 0.5*1e-8*x
u = 3*1e-9*x
g = 0*1e-8*x
k = 440
h = 0.5
f0 = 0.1
f1 = 0.2
f2 = 0.3
f3 = 0.4
gamma_cutoff = 5
t0 = 0.0
t1 = h*(1/(2*Nanc))
t1half = h*(gamma_cutoff/(2*Nanc))
t2 = h*(10/(2*Nanc))
t3 = h*(100/(2*Nanc))
t4 = h*1.0
mean, shape, proportion_synonymous = 500, 0.5, 0.3
""")
        params_path = params_file.name
    
    try:
        # Create output file path
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as out_file:
            out_path = out_file.name
        
        # Run bvalcalc with --genome and extra flags
        cmd = [
            "python", "-m", "Bvalcalc", 
            "--genome",
            "--params", params_path,
            "--bedgff", bed_path,
            "--rec_map", rec_path,
            "--chr_sizes", chr_path,
            "--out", out_path,
            "--out_binsize", "1000",
            "--quiet"  # Keep output clean
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Check that command succeeded
        assert result.returncode == 0, f"Command failed: {result.stderr}"
        
        # Check that output file was created
        assert os.path.exists(out_path), "Output file was not created"
        
        # Read the output file
        with open(out_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        # Verify we have at least the header lines
        assert len(lines) >= 4, f"Expected at least 4 lines (header + data), got {len(lines)}"
        
        # Test first line: version (should match current package version)
        from Bvalcalc import __version__
        expected_version = f"# Bvalcalc v{__version__}"
        assert lines[0] == expected_version, f"Expected version line '{expected_version}', got: {lines[0]}"
        
        # Test second line: command (should contain our flags)
        command_line = lines[1]
        assert command_line.startswith("# "), f"Command line should start with '# ', got: {command_line}"
        assert "--genome" in command_line, "Command line should contain --genome"
        assert "--params" in command_line, "Command line should contain --params"
        assert "--bedgff" in command_line, "Command line should contain --bedgff"
        assert "--rec_map" in command_line, "Command line should contain --rec_map"
        assert "--chr_sizes" in command_line, "Command line should contain --chr_sizes"
        assert "--out_binsize" in command_line, "Command line should contain --out_binsize"
        
        # Check if we have a warning line (may or may not be present depending on rec_map coverage)
        if lines[2].startswith("# WARNING:"):
            # Warning is present - check it's about start coverage
            warning_line = lines[2]
            assert "doesn't cover start of chromosome" in warning_line, "Warning should mention start coverage"
            # Format should be on line 3
            format_line = lines[3]
            assert format_line == "# Format: Chromosome,Start,B", f"Expected format line, got: {format_line}"
            data_start_idx = 4
        else:
            # No warning - format should be on line 2
            format_line = lines[2]
            assert format_line == "# Format: Chromosome,Start,B", f"Expected format line, got: {format_line}"
            data_start_idx = 3
        
        # Test that data starts after header
        assert len(lines) > data_start_idx, f"Should have data lines after header (starting at index {data_start_idx})"
        data_line = lines[data_start_idx]
        assert data_line.startswith("chr1,"), f"Data line should start with 'chr1,', got: {data_line}"
        
        print("✓ Header format test passed!")
        print(f"  Version: {lines[0]}")
        print(f"  Command: {lines[1][:50]}...")
        if data_start_idx == 4:
            print(f"  Warning: {lines[2][:50]}...")
            print(f"  Format: {lines[3]}")
            print(f"  Data: {lines[4]}")
        else:
            print(f"  Format: {lines[2]}")
            print(f"  Data: {lines[3]}")
        
    finally:
        # Clean up temporary files
        for path in [bed_path, rec_path, chr_path, params_path, out_path]:
            if os.path.exists(path):
                os.unlink(path)


if __name__ == "__main__":
    test_header_output_format()
