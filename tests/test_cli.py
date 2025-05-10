import subprocess
import sys
import re
from pathlib import Path

# Invoke the CLI via the installed package entry point.
BASE_CMD = [sys.executable, "-m", "bvalcalc.cli"]

def test_cli_site_basic():
    # python -m bvalcalc.cli --site --pop_params tests/testparams/gcBasicParams.py --distance 100 --gene_size 5000
    params = Path(__file__).parent / "testparams" / "gcBasicParams.py"
    cmd = BASE_CMD + [
        "--site",
        "--pop_params", str(params),
        "--distance", "100",
        "--gene_size", "5000",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "B for site 100bp away from 5000bp region:" in result.stdout
    assert "0.9246145075198539" in result.stdout

def test_cli_gene_basic():
    # python -m bvalcalc.cli --gene --pop_params tests/testparams/nogcBasicParams.py --gene_size 10000
    params = Path(__file__).parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    cmd = BASE_CMD + [
        "--gene",
        "--pop_params", str(params),
        "--gene_size", "10000",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "====== R E S U L T S ! =============================" in out
    assert "B for adjacent site: 0.8049242606161049" in out
    assert "Mean B for flanking region: 0.9761402805820517" in out
    assert "No output CSV requested; skipping save." in out
    assert "= B value calculated" in out

def test_cli_gene_gcparams(tmp_path):
    # python -m bvalcalc.cli --gene --pop_params tests/testparams/gcBasicParams.py --gene_size 10000 --plot_output tests/testout/test_plot.png
    params   = Path(__file__).parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    plot_path = Path(__file__).parents[1] / "tests" / "testout" / "test_plot.png"
    cmd = BASE_CMD + [
        "--gene",
        "--pop_params", str(params),
        "--gene_size", "10000",
        "--plot_output", str(plot_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "B for adjacent site: 0.8910346781386976" in out
    assert "Mean B for flanking region: 0.9810661565709757" in out
    assert "B at start and end of the neutral region" in out
    assert "====== P L O T T I N G . . . =======================" in out
    assert f"Plot saved to {plot_path}" in out
    assert plot_path.exists(), f"Expected plot at {plot_path}, but not found"

def test_cli_genome_basic(tmp_path):
    # python -m bvalcalc.cli --genome --pop_params tests/testparams/nogcBasicParams.py --bedgff_path tests/testfiles/200kb_slimtest.csv --chr_sizes examples/test_sizes.txt --out <tmp>/200kb_dfe5.bvals --out_binsize 1
    params         = Path(__file__).parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    chr_sizes_path = Path(__file__).parents[1] / "examples" / "test_sizes.txt"
    output_path    = tmp_path / "200kb_dfe5.bvals"
    cmd = BASE_CMD + [
        "--genome",
        "--pop_params", str(params),
        "--bedgff_path", str(bed_path),
        "--chr_sizes", str(chr_sizes_path),
        "--out", str(output_path),
        "--out_binsize", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.753693843332109" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_genome_gcparams(tmp_path):
    # python -m bvalcalc.cli --genome --pop_params tests/testparams/gcBasicParams.py --bedgff_path tests/testfiles/200kb_slimtest.csv --chr_sizes examples/test_sizes.txt --out <tmp>/gc_bvals.bvals --out_binsize 1
    params         = Path(__file__).parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    bed_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    chr_sizes_path = Path(__file__).parents[1] / "examples" / "test_sizes.txt"
    output_path    = tmp_path / "gc_bvals.bvals"
    cmd = BASE_CMD + [
        "--genome",
        "--pop_params", str(params),
        "--bedgff_path", str(bed_path),
        "--chr_sizes", str(chr_sizes_path),
        "--out", str(output_path),
        "--out_binsize", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "====== R E S U L T S ====== S U M M A R Y ==========" in out
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.836347850423207" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_genome_with_recmap_plot(tmp_path):
    # python -m bvalcalc.cli --genome --pop_params tests/testparams/nogcBasicParams.py --bedgff_path tests/testfiles/200kb_slimtest.csv --chr_sizes examples/test_sizes.txt --rec_map tests/testfiles/200kb.map --out <tmp>/200kb_dfe5.bvals --out_binsize 1
    params         = Path(__file__).parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    map_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb.map"
    chr_sizes_path = Path(__file__).parents[1] / "examples" / "test_sizes.txt"
    output_path    = tmp_path / "200kb_dfe5.bvals"
    cmd = BASE_CMD + [
        "--genome",
        "--pop_params", str(params),
        "--bedgff_path", str(bed_path),
        "--chr_sizes", str(chr_sizes_path),
        "--rec_map", str(map_path),
        "--out", str(output_path),
        "--out_binsize", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Cumulative length of chromosome under selection: 99990bp (50.0%)" in out
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.7015847245703709" in out
    assert f"Appended B values to: {output_path.as_posix()}" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_mean_b_value(tmp_path):
    # python -m bvalcalc.cli --region --pop_params tests/testparams/nogcBasicParams.py --bedgff_path tests/testfiles/200kb_slimtest.csv --plot_output tests/testout/genome_test.png --calc_region chr_200kb:1514-62456
    params = "tests/testparams/nogcBasicParams.py"
    cmd = BASE_CMD + [
        "--region",
        "--pop_params", params,
        "--bedgff_path", "tests/testfiles/200kb_slimtest.csv",
        "--plot_output", "tests/testout/genome_test.png",
        "--calc_region", "chr_200kb:1514-62456",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"Process failed: {result.stderr}"
    assert "Plot saved to tests/testout/genome_test.png" in out
    match = re.search(r"Mean B of neutral sites across specified region:\s+([0-9.]+)", result.stdout)
    assert match, "Could not find mean B output in CLI output"
    mean_b = float(match.group(1))
    expected = 0.7609515711751818
    assert abs(mean_b - expected) < 1e-10, f"Expected {expected}, got {mean_b}"

def test_cli_gene_contract():
    # python -m bvalcalc.cli --gene --pop_params tests/testparams/ContractParams_5N_1T.py --pop_change --plot_output
    params = Path(__file__).parents[1] / "tests" / "testparams" / "ContractParams_5N_1T.py"
    cmd = BASE_CMD + [
        "--gene",
        "--pop_params", str(params),
        "--pop_change",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B for flanking region: 0.9623874578845304" in out

def test_cli_gene_expand():
    # python -m bvalcalc.cli --gene --pop_params tests/testparams/ExpandParams_5N_1T.py --pop_change
    params = Path(__file__).parents[1] / "tests" / "testparams" / "ExpandParams_5N_1T.py"
    cmd = BASE_CMD + [
        "--gene",
        "--pop_params", str(params),
        "--pop_change",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B for flanking region: 0.9823767432235562" in out
    assert "B prior to demographic change" in out
    assert "B post B-calculation" in out

def test_cli_selfing():
    # python -m bvalcalc.cli --gene --pop_params tests/testparams/SelfParams_0.9S_0.5h.py
    params = Path(__file__).parents[1] / "tests" / "testparams" / "SelfParams_0.9S_0.5h.py"
    cmd = BASE_CMD + [
        "--gene",
        "--pop_params", str(params),
        "--pop_change",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "B for adjacent site: 0.6034770660828896" in out
    assert "Mean B for flanking region: 0.8043714716235398" in out
    assert "B at start and end of the neutral region: [0.60347707 0.6035028  0.60352854 ... 0.89001725 0.89001929 0.89002133]" in out



