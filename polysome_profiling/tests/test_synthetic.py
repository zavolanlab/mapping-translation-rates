"""
Synthetic-data test for polysome_profiling/polysome_analysis_model_RL_pool.py.

synth.py (same directory) generates a count matrix from the model stated in the Methods
(ribosomes per molecule ~ Poisson(λ_i); the same amount of spike-in in every
fraction; reads ∝ molecules × length × sequencing depth of the fraction). The
true λ_i are therefore known by construction, and they are the oracle: whatever
implements the Methods must recover them, whatever the other genes look like
and however deeply each fraction was sequenced.

The tests write the three input files into a temporary directory, run the
script there, and read its output. They need no data files and take about 7 s:
    python -m pytest polysome_profiling/tests
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import synth  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "polysome_analysis_model_RL_pool.py"
TRUE_LAMBDA = np.array([0.5, 1.0, 2.0, 3.5, 5.0, 7.0])


def run_script(tmp_path, counts, len_df):
    """Write the three inputs the script reads, run it, return its protein-coding output."""
    counts.to_csv(tmp_path / "all_rna_on_exon_reads.tsv", sep="\t", index=False)
    with open(tmp_path / "transcriptome.fa", "w") as fa:
        for _, r in len_df.iterrows():
            tag = " CDS=1-30" if r.protein_coding else ""
            fa.write(f">{r.transcript_ID}{tag}\n" + "A" * int(r.transcript_length) + "\n")
    len_df[["gene_ID", "transcript_ID"]].to_csv(tmp_path / "gene_id_transcript_id.tsv", sep="\t", index=False, header=False)
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT)], cwd=tmp_path, check=True, capture_output=True)
    return pd.read_csv(tmp_path / "pc_frequencies_mrl_poisson_fit.tsv", sep="\t").set_index("gene_ID")


def test_recovers_known_lambda(tmp_path):
    """Six genes with known λ in a background of 300 genes of varied λ (so that the
    mRNA content differs between fractions, as in a real gradient). Noise-free counts:
    the recovered λ must equal the true λ to better than 2 %."""
    rng = np.random.default_rng(1)
    background = rng.uniform(0.3, 9.0, 300)
    lam = np.concatenate([TRUE_LAMBDA, background])
    counts, len_df = synth.make_dataset(lam, n_molecules=np.full(len(lam), 1e4))
    out = run_script(tmp_path, counts, len_df)
    got = out.loc[len_df.gene_ID[: len(TRUE_LAMBDA)], "ribosome_load"].values
    rel = np.abs(got - TRUE_LAMBDA) / TRUE_LAMBDA
    assert rel.max() < 0.02, f"true λ {TRUE_LAMBDA}, recovered {np.round(got, 2)}"


def test_invariant_to_sequencing_depth_of_each_fraction(tmp_path):
    """Sequencing some fractions 10× deeper than others must not change λ: the spike-in
    absorbs the depth (that is what it is for)."""
    rng = np.random.default_rng(2)
    lam = np.concatenate([TRUE_LAMBDA, rng.uniform(0.3, 9.0, 100)])
    samples = [*synth.FRACTIONS, synth.HEAVY]
    depth = {s: d for s, d in zip(samples, rng.uniform(0.2, 2.0, len(samples)))}
    c1, l1 = synth.make_dataset(lam)
    c2, l2 = synth.make_dataset(lam, depth=depth)
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    a = run_script(tmp_path / "a", c1, l1)
    b = run_script(tmp_path / "b", c2, l2)
    ids = l1.gene_ID[: len(TRUE_LAMBDA)]
    assert np.allclose(a.loc[ids, "ribosome_load"], b.loc[ids, "ribosome_load"], rtol=1e-3), (
        f"λ changed with sequencing depth: {np.round(a.loc[ids, 'ribosome_load'].values, 2)} vs "
        f"{np.round(b.loc[ids, 'ribosome_load'].values, 2)}")


def test_column_order_is_irrelevant(tmp_path):
    """Metamorphic test: shuffling the columns of the count matrix must not change any
    result, since the Methods identify fractions by name, not by position."""
    rng = np.random.default_rng(3)
    lam = np.concatenate([TRUE_LAMBDA, rng.uniform(0.3, 9.0, 100)])
    counts, len_df = synth.make_dataset(lam)
    cols = list(counts.columns[1:]); rng.shuffle(cols)
    shuffled = counts[["gene_ID", *cols]]
    (tmp_path / "a").mkdir(); (tmp_path / "b").mkdir()
    a = run_script(tmp_path / "a", counts, len_df)
    b = run_script(tmp_path / "b", shuffled, len_df)
    ids = len_df.gene_ID[: len(TRUE_LAMBDA)]
    assert np.allclose(a.loc[ids, "ribosome_load"], b.loc[ids, "ribosome_load"]), "λ depends on the column order of the count matrix"
