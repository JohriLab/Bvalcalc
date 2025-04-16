import subprocess
import sys
import re
from pathlib import Path

def test_cli_site_basic():
    #python Bvalcalc.py --site --pop_params tests/testparams/gcBasicParams.py --distance 100 --gene_size 5000
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

def test_cli_gene_basic():
    # ./Bvalcalc.py --gene --pop_params tests/testparams/nogcBasicParams.py --gene_size 10000
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "nogcBasicParams.py"

    result = subprocess.run(
        [sys.executable, str(script),
         "--gene",
         "--pop_params", str(params),
         "--gene_size", "10000"],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    out = result.stdout + result.stderr

    assert "====== R E S U L T S ! =============================" in out
    assert "B for adjacent site: 0.8049242606161049" in out
    assert "Mean B for flanking region: 0.9761402805820517" in out
    assert "No output CSV requested; skipping save." in out
    assert "= B value calculated" in out


def test_cli_gene_gcparams(tmp_path):
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    output_path = tmp_path / "40kb_gc.bvals"

    result = subprocess.run(
        [sys.executable, str(script),
         "--gene",
         "--pop_params", str(params),
         "--gene_size", "10000",
         "--out", str(output_path),
         "--plot_output"],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    
    out = result.stdout + result.stderr
    assert "= Calculating relative diversity (B) for a neutral region adjacent to a single selected element =" in out
    assert "Length of element under selection: 10000bp" in out
    assert "Length of flanking neutral region: 40000bp" in out
    assert "====== S T A R T I N G ===== C A L C ===============" in out
    assert "====== F I N I S H E D ===== C A L C ===============" in out
    assert "====== R E S U L T S ! =============================" in out
    assert "B for adjacent site: 0.8910346781386976" in out
    assert "Mean B for flanking region: 0.9810661565709757" in out
    assert "B at start and end of the neutral region" in out
    assert "====== P L O T T I N G . . . =======================" in out
    assert "Plot saved to Bplot.png" in out
    assert f"Saved B values to: {output_path.as_posix()}" in out
    assert "= B value calculated" in out

def test_cli_genome_basic(tmp_path):
    #python Bvalcalc.py --genome --pop_params tests/testparams/nogcBasicParams.py --bedgff_path exampleData/200kb_slimtest.csv --chr_size 200000 --out /path/to/output/dgas.bvals
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path = Path(__file__).resolve().parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    chr_sizes_path = Path(__file__).resolve().parents[1] / "exampleData" / "test_sizes.txt"
    output_path = tmp_path / "200kb_dfe5.bvals"

    result = subprocess.run(
        [sys.executable, str(script),
            "--genome",
            "--pop_params", str(params),
            "--bedgff_path", str(bed_path),
            "--chr_size", str(chr_sizes_path),
            "--out", str(output_path)],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    out = result.stdout + result.stderr
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.753693843332109" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_genome_gcparams(tmp_path):
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    bed_path = Path(__file__).resolve().parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    chr_sizes_path = Path(__file__).resolve().parents[1] / "exampleData" / "test_sizes.txt"
    output_path = tmp_path / "gc_bvals.bvals"

    result = subprocess.run(
        [sys.executable, str(script),
         "--genome",
         "--pop_params", str(params),
         "--bedgff_path", str(bed_path),
         "--chr_sizes", str(chr_sizes_path),
         "--out", str(output_path)],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    out = result.stdout + result.stderr
    assert "====== R E S U L T S ====== S U M M A R Y ==========" in out
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.836347850423207" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_genome_with_recmap_plot(tmp_path):
    script = Path(__file__).resolve().parents[1] / "Bvalcalc.py"
    params = Path(__file__).resolve().parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path = Path(__file__).resolve().parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    map_path = Path(__file__).resolve().parents[1] / "tests" / "testfiles" / "200kb.map"
    chr_sizes_path = Path(__file__).resolve().parents[1] / "exampleData" / "test_sizes.txt"
    output_path = tmp_path / "200kb_dfe5.bvals"

    result = subprocess.run(
        [sys.executable, str(script),
         "--genome",
         "--pop_params", str(params),
         "--bedgff_path", str(bed_path),
         "--chr_sizes", str(chr_sizes_path),
         "--rec_map", str(map_path),
         "--out", str(output_path)],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    out = result.stdout + result.stderr
    assert "Cumulative length of chromosome under selection: 99990bp (50.0%)" in out
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.7015847245703709" in out
    assert f"Saved B values to: {output_path.as_posix()}" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_mean_b_value():
    #./Bvalcalc.py --region --pop_params tests/testparams/nogcBasicParams.py --bedgff_path exampleData/200kb_slimtest.csv --chr_start 1 --chr_size 200000 --plot_output --calc_start 1514 --calc_end 62456

    result = subprocess.run([
        sys.executable, "Bvalcalc.py",
        "--region",
        "--pop_params", "tests/testparams/nogcBasicParams.py",
        "--bedgff_path", "tests/testfiles/200kb_slimtest.csv",
        "--plot_output",
        "--calc_region", "chr_200kb:1514-62456",
    ], capture_output=True, text=True)

    assert result.returncode == 0, f"Process failed: {result.stderr}"

    match = re.search(r"Mean B of neutral sites across specified region:\s+([0-9.]+)", result.stdout)
    assert match, "Could not find mean B output in CLI output"

    mean_b = float(match.group(1))
    expected = 0.7609515711751818
    assert abs(mean_b - expected) < 1e-10, f"Expected {expected}, got {mean_b}"