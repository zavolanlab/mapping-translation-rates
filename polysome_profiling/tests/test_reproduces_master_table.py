"""
Regression test for polysome_profiling/polysome_analysis_model_RL_pool.py:
the script must regenerate the mean ribosome loads published with
González et al., Sci Data (2025), Zenodo record 10.5281/zenodo.15858175.

Place it at  polysome_profiling/tests/test_reproduces_master_table.py  and run
    MRL_DATA_DIR=/path/with/inputs  python -m pytest polysome_profiling/tests
where MRL_DATA_DIR contains the three inputs the script reads
(all_rna_on_exon_reads.tsv, transcriptome.fa incl. the spike-in sequences,
gene_id_transcript_id.tsv incl. the spike-ins) and the master table exported as
published_mrl.tsv (columns gene_ID, ribosome_load, error_ribosome_load).
The test is skipped when MRL_DATA_DIR is not set, so the repository's CI stays green
without the 80 MB of inputs; run it locally before every change to the script.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "polysome_profiling" / "polysome_analysis_model_RL_pool.py"
DATA = os.environ.get("MRL_DATA_DIR")


@pytest.fixture(scope="module")
def regenerated(tmp_path_factory):
    if not DATA:
        pytest.skip("set MRL_DATA_DIR to the directory holding the script inputs and published_mrl.tsv")
    needed = ("all_rna_on_exon_reads.tsv", "transcriptome.fa", "gene_id_transcript_id.tsv", "published_mrl.tsv")
    missing = [n for n in needed if not (Path(DATA) / n).exists()]
    if missing:
        pytest.skip(f"MRL_DATA_DIR={DATA} lacks {', '.join(missing)}")
    work = tmp_path_factory.mktemp("mrl")
    for name in needed[:3]:
        shutil.copy(Path(DATA) / name, work / name)
    subprocess.run([sys.executable, str(SCRIPT)], cwd=work, check=True)
    ours = pd.read_csv(work / "pc_frequencies_mrl_poisson_fit.tsv", sep="\t").set_index("gene_ID")
    pub = pd.read_csv(Path(DATA) / "published_mrl.tsv", sep="\t").set_index("gene_ID")
    return pub.join(ours, lsuffix="_pub", rsuffix="_new", how="inner"), set(pub.index) - set(ours.index)


def test_all_published_genes_are_produced(regenerated):
    _, missing = regenerated
    assert not missing, f"{len(missing)} published genes missing from the output"


@pytest.mark.parametrize("col", ["ribosome_load", "error_ribosome_load"])
def test_column_reproduces(regenerated, col):
    both, _ = regenerated
    a, b = both[f"{col}_pub"], both[f"{col}_new"]
    ok = np.isfinite(a) & np.isfinite(b)
    r = np.corrcoef(a[ok], b[ok])[0, 1]
    assert r > 0.99, f"{col}: Pearson r = {r:.4f} between published and regenerated values ({ok.sum()} genes)"
    if col == "ribosome_load":
        med = float(np.median(np.abs(a - b)))
        assert med < 0.05, f"median |difference| in mean ribosome load = {med:.3f}"
