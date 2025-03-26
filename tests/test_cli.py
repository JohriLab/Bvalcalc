import subprocess
import sys
from pathlib import Path

def test_cli_site_basic():
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parent / "testparams" / "LowRParams.py"

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
    assert "0.9294346764091733" in result.stdout
