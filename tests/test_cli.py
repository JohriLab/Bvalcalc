import subprocess
import sys
import re
from pathlib import Path

# Invoke the CLI via the installed package entry point.
BASE_CMD = [sys.executable, "-m", "Bvalcalc.cli"]

def test_cli_site_basic():
    # poetry run Bvalcalc --site --params tests/testparams/gcBasicParams.py --distance 100 --element_size 5000
    params = Path(__file__).parent / "testparams" / "gcBasicParams.py"
    cmd = BASE_CMD + [
        "--site",
        "--params", str(params),
        "--distance", "100",
        "--element_size", "5000",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "B for site 100bp away from 5000bp region:" in result.stdout
    assert "0.9246145075198539" in result.stdout

def test_cli_gene_basic():
    # poetry run Bvalcalc --gene --params tests/testparams/nogcBasicParams.py --element_size 10000
    params = Path(__file__).parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    cmd = BASE_CMD + [
        "--gene",
        "--params", str(params),
        "--element_size", "10000",
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
    # poetry run Bvalcalc --gene --params tests/testparams/gcBasicParams.py --element_size 10000 --plot tests/testout/test_plot.png
    params   = Path(__file__).parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    plot_path = Path(__file__).parents[1] / "tests" / "testout" / "test_plot.png"
    cmd = BASE_CMD + [
        "--gene",
        "--params", str(params),
        "--element_size", "10000",
        "--plot", str(plot_path),
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
    # poetry run Bvalcalc --genome --params tests/testparams/nogcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --chr_sizes tests/testfiles/test_sizes.txt
    # poetry run Bvalcalc --region chr_200kb:1-200000 --params tests/testparams/nogcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --plot
    params         = Path(__file__).parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    chr_sizes_path = Path(__file__).parents[1] / "tests" / "testfiles" / "test_sizes.txt"
    output_path    = tmp_path / "200kb_dfe5.bvals"
    cmd = BASE_CMD + [
        "--genome",
        "--params", str(params),
        "--bedgff", str(bed_path),
        "--chr_sizes", str(chr_sizes_path),
        "--out", str(output_path),
        "--out_binsize", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    match = re.search(r"Mean B of neutral sites across chromosome chr_200kb: ([0-9.]+)", out)
    assert match, "Could not find mean B output in CLI output"
    mean_b = float(match.group(1))
    expected = 0.753693843332108
    assert abs(mean_b - expected) < 1e-10, f"Expected {expected}, got {mean_b}"
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_genome_gcparams(tmp_path):
    # poetry run Bvalcalc --genome --params tests/testparams/gcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --chr_sizes tests/testfiles/test_sizes.txt --out <tmp>/gc_bvals.bvals --out_binsize 1
    params         = Path(__file__).parents[1] / "tests" / "testparams" / "gcBasicParams.py"
    bed_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    chr_sizes_path = Path(__file__).parents[1] / "tests" / "testfiles" / "test_sizes.txt"
    output_path    = tmp_path / "gc_bvals.bvals"
    cmd = BASE_CMD + [
        "--genome",
        "--params", str(params),
        "--bedgff", str(bed_path),
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
    # poetry run Bvalcalc --genome --params tests/testparams/nogcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --chr_sizes tests/testfiles/test_sizes.txt --rec_map tests/testfiles/200kb.map --out examples/200kb_dfe5.bvals --out_binsize 1
    params         = Path(__file__).parents[1] / "tests" / "testparams" / "nogcBasicParams.py"
    bed_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_slimtest.csv"
    map_path       = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb.map"
    chr_sizes_path = Path(__file__).parents[1] / "tests" / "testfiles" / "test_sizes.txt"
    output_path    = tmp_path / "200kb_dfe5.bvals"
    cmd = BASE_CMD + [
        "--genome",
        "--params", str(params),
        "--bedgff", str(bed_path),
        "--chr_sizes", str(chr_sizes_path),
        "--rec_map", str(map_path),
        "--out", str(output_path),
        "--out_binsize", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Cumulative length of chromosome under selection: 99990bp (50.0%)" in out
    match = re.search(r"Mean B of neutral sites across chromosome chr_200kb: ([0-9.]+)", out)
    assert match, "Could not find mean B output in CLI output"
    mean_b = float(match.group(1))
    expected = 0.701584724570370
    assert abs(mean_b - expected) < 1e-10, f"Expected {expected}, got {mean_b}"
    assert f"Appended B values to: {output_path.as_posix()}" in out
    assert output_path.exists(), "Expected output file not created"
    assert output_path.stat().st_size > 0, "Output file is empty"

def test_cli_mean_b_value(tmp_path):
    # poetry run Bvalcalc --region chr_200kb:1514-62456 --params tests/testparams/nogcBasicParams.py --bedgff tests/testfiles/200kb_slimtest.csv --plot tests/testout/genome_test.png 
    params = "tests/testparams/nogcBasicParams.py"
    cmd = BASE_CMD + [
        "--region", "chr_200kb:1514-62456",
        "--params", params,
        "--bedgff", "tests/testfiles/200kb_slimtest.csv",
        "--plot", "tests/testout/genome_test.png",
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
    # poetry run Bvalcalc --gene --params tests/testparams/ContractParams_5N_1T.py --pop_change --plot
    params = Path(__file__).parents[1] / "tests" / "testparams" / "ContractParams_5N_1T.py"
    cmd = BASE_CMD + [
        "--gene",
        "--params", str(params),
        "--pop_change",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B for flanking region: 0.9623874578845304" in out

def test_cli_gene_expand():
    # poetry run Bvalcalc --gene --params tests/testparams/ExpandParams_5N_1T.py --pop_change
    params = Path(__file__).parents[1] / "tests" / "testparams" / "ExpandParams_5N_1T.py"
    cmd = BASE_CMD + [
        "--gene",
        "--params", str(params),
        "--pop_change",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B for flanking region: 0.982376743223556" in out
    assert "B prior to demographic change" in out
    assert "B post B-calculation" in out

def test_cli_selfing():
    # poetry run Bvalcalc --gene --params tests/testparams/SelfParams_0.9S_0.5h.py
    params = Path(__file__).parents[1] / "tests" / "testparams" / "SelfParams_0.9S_0.5h.py"
    cmd = BASE_CMD + [
        "--gene",
        "--params", str(params),
        "--pop_change",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "B for adjacent site: 0.6034770660828896" in out
    assert "Mean B for flanking region: 0.8043714716235398" in out
    assert "B at start and end of the neutral region: [0.60347707 0.6035028  0.60352854 ... 0.89001725 0.89001929 0.89002133]" in out

def test_cli_positions_minimum_filter():
    # poetry run Bvalcalc --Bmap ./tests/testfiles/false_Bvalues_chr3R.csv --positions ./tests/testfiles/posfile.csv --out_minimum 0.5
    bmap_path = Path(__file__).parents[1] / "tests" / "testfiles" / "false_Bvalues_chr3R.csv"
    pos_path  = Path(__file__).parents[1] / "tests" / "testfiles" / "posfile.csv"

    cmd = BASE_CMD + [
        "--Bmap", str(bmap_path),
        "--positions", str(pos_path),
        "--out_minimum", "0.5",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B across filtered sites: 0.562500" in out
    assert "Max B across filtered sites: 0.750000 at chr_2R:20000" in out
    assert "Min B across filtered sites: 0.500000 at chr_2R:1" in out

def test_cli_genome_hri_marking(tmp_path):
    # poetry run Bvalcalc --region chr_200kb:1-199999  --params ./tests/testparams/InterfParams.py  --bedgff ./tests/testfiles/200kb_hri.csv --rec_map ./tests/testfiles/200kb_hri.map --precise_chunks 2 --out haut.B --out_binsize 100 --hri --plot
    params      = Path(__file__).parents[1] / "tests" / "testparams" / "InterfParams.py"
    bed_path    = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_hri.csv"
    map_path    = Path(__file__).parents[1] / "tests" / "testfiles" / "200kb_hri.map"
    output_path = Path(__file__).parents[1] / "tests" / "testout" / "hri_200kb_test.B"

    cmd = BASE_CMD + [
        "--genome",
        "--params", str(params),
        "--bedgff", str(bed_path),
        "--rec_map", str(map_path),
        "--precise_chunks", "2",
        "--out", str(output_path),
        "--out_binsize", "100",
        "--hri",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    out = result.stdout + result.stderr

    # CLI ran
    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert "Mean B of neutral sites across chromosome chr_200kb: 0.588623210410594" in out

    # File wrote
    assert output_path.exists(), "Expected output file not created"

    # Contents
    with open(output_path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    assert len(lines) == 2003, f"Expected 2003 lines (3 header + 2000 data), found {len(lines)}"

    # Specific rows must be present
    assert "chr_200kb,1,0.340332'" in lines
    assert "chr_200kb,101,0.340332'" in lines
    assert "chr_200kb,40101,0.359374'" in lines

    # Final line check
    assert lines[-1] == "chr_200kb,199901,0.327636'", f"Unexpected final line: {lines[-1]}"

