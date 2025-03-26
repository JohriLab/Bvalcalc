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
