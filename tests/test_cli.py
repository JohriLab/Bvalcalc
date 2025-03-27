import subprocess
import sys
from pathlib import Path

def test_cli_site_basic():
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parent / "testparams" / "gcBasicParams.py"

    result = subprocess.run(
        [sys.executable, str(script),
         "--site",
         "--pop_params", str(params),
         "--distance", "100",
         "--gene_size", "5000"],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "B for site 100bp away from 5000bp region:" in result.stdout
    assert "0.9246145075198539" in result.stdout

def test_cli_region_basic(tmp_path):
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    output_path = tmp_path / "dfe_bvals.csv"  # Use tmp_path to avoid clobbering real files

    result = subprocess.run(
        [sys.executable, str(script),
         "--region",
         "--pop_params", str(params),
         "--gene_size", "10000",
         "--out", str(output_path),
         "--silent"],
        capture_output=True,
        text=True
    )

    # Check CLI executed correctly
    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"

    # Check output file exists and is not empty
    assert output_path.exists(), "Expected output CSV was not created"
    assert output_path.stat().st_size > 0, "Output CSV is empty"

    # Optional: check specific content
    contents = output_path.read_text()
    assert "Distance,B" in contents, "Output CSV does not contain expected header"

def test_cli_region_gcparams(tmp_path):
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    output_path = tmp_path / "40kb_gc.bvals"

    result = subprocess.run(
        [sys.executable, str(script),
         "--region",
         "--pop_params", str(params),
         "--gene_size", "10000",
         "--out", str(output_path),
         "--plot_output"],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    
    out = result.stdout
    assert "= Calculating relative diversity (B) for a neutral region adjacent to a single selected element =" in out
    assert "Distribution of fitness effects (DFE): 40000bp" in out
    assert "Length of region under selection: 10000bp" in out
    assert "Length of flanking neutral region: 40000bp" in out
    assert "====== S T A R T I N G ===== C A L C ===============" in out
    assert "====== F I N I S H E D ===== C A L C ===============" in out
    assert "====== R E S U L T S ! =============================" in out
    assert "B for adjacent site: 0.8910346781386976" in out
    assert "Mean B for flanking region: 0.9810661565709757" in out
    assert "B at start and end of the neutral region" in out
    assert "====== P L O T T I N G . . . =======================" in out
    assert "Plot saved to region_plot.png" in out
    assert f"Saved B values to: {output_path}" in out
    assert "= B value calculated" in out

def test_cli_genome_basic(tmp_path):
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path = Path(__file__).resolve().parents[1] / "exampleData" / "200kb_slimtest.csv"
    output_path = tmp_path / "200kb_dfe5.bvals"

    result = subprocess.run(
        [sys.executable, str(script),
            "--genome",
            "--pop_params", str(params),
            "--bedgff_path", str(bed_path),
            "--chr_start", "1",
            "--chr_end", "200000",
            "--out", str(output_path)],  # <- removed --silent
        capture_output=True,
        text=True
    )


    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"

    out = result.stdout
    assert "Cumulative length of regions under selection: 106994bp" in out
    assert "Mean B of neutral sites across genome: 0.738525356400387" in out
