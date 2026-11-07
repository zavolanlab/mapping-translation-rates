"""
Synthetic polysome-profiling count matrices in the exact input format the
authors' script expects, generated from the model stated in the paper's Methods:

  * each mRNA i has a mean ribosome load lambda_i; the number of ribosomes k on a
    molecule is Poisson(k; lambda_i)
  * fractions resolve k = 0..10 ribosomes (free/40S -> 0, 80S -> 1, disome -> 2, ...),
    everything with k > 10 ends up in the 'heavysome' fraction
  * the same amount S of spike-in RNA is added to every fraction (so spike-in
    molecules are constant across fractions)
  * reads for gene i in fraction k are proportional to molecules x transcript
    length x fraction-specific sequencing depth

Because the generator implements the *Methods*, any analysis that implements
the Methods must recover lambda_i. That is the oracle.
"""
import numpy as np
import pandas as pd
from scipy.stats import poisson

# Sample names exactly as hard-coded in the authors' sample_dict / heavy_fract
FRACTIONS = {
    "CDS_RNAseq_Sample_HepG2polysomeFreeFraction": 0,
    "CDS_RNAseq_Sample_HepG2polysome40s": 0,
    "CDS_RNAseq_Sample_HepG2polysome80s": 1,
    "CDS_RNAseq_Sample_HepG2polysomedisome": 2,
    "CDS_RNAseq_Sample_HepG2polysometrisome": 3,
    "CDS_RNAseq_Sample_HepG2polysometetrasome": 4,
    "CDS_RNAseq_Sample_HepG2polysomepentasome": 5,
    "CDS_RNAseq_Sample_HepG2polysomehexasome": 6,
    "CDS_RNAseq_Sample_HepG2polysomeheptasome": 7,
    "CDS_RNAseq_Sample_HepG2polysomeoctasome": 8,
    "CDS_RNAseq_Sample_HepG2polysomenonasome": 9,
    "CDS_RNAseq_Sample_HepG2polysomedecasome": 10,
}
HEAVY = "CDS_RNAseq_Sample_HepG2polysomeheavysome"
KMAX = 10


def molecules_per_fraction(lam, n_molecules):
    """Expected molecules of one mRNA species in each sample (dict name->molecules)."""
    pk = poisson.pmf(np.arange(KMAX + 1), lam)
    tail = 1.0 - pk.sum()
    out = {}
    for name, k in FRACTIONS.items():
        # free fraction and 40S both carry k=0 molecules; split 50/50 between them
        share = 0.5 if k == 0 else 1.0
        out[name] = n_molecules * pk[k] * share
    out[HEAVY] = n_molecules * tail
    return out


def make_dataset(
    lambdas,
    lengths=None,
    n_molecules=None,
    spikein_molecules=1e5,
    spikein_length=1000,
    depth=None,
    replicates=("rep1", "rep2", "rep3"),
    gene_prefix="ENSG",
    noise=False,
    seed=0,
    reads_per_molecule_nt=1e-3,
):
    """
    Returns (count_matrix, len_df) in the authors' input format.

    lambdas:  array of true mean ribosome loads, one per gene
    lengths:  transcript lengths (nt); default 2000
    n_molecules: molecules per gene in the lysate; default 1e4
    depth:    dict sample -> sequencing-depth multiplier (default 1 everywhere)
    noise:    if True, Poisson-sample the read counts
    """
    rng = np.random.default_rng(seed)
    lambdas = np.asarray(lambdas, dtype=float)
    n = len(lambdas)
    lengths = np.full(n, 2000.0) if lengths is None else np.asarray(lengths, float)
    n_molecules = np.full(n, 1e4) if n_molecules is None else np.asarray(n_molecules, float)
    samples = [*FRACTIONS.keys(), HEAVY]
    depth = {s: 1.0 for s in samples} if depth is None else depth

    gene_ids = [f"{gene_prefix}{i:011d}" for i in range(n)]
    rows = {}
    for gid, lam, L, N in zip(gene_ids, lambdas, lengths, n_molecules):
        mols = molecules_per_fraction(lam, N)
        rows[gid] = {s: mols[s] * L * depth[s] * reads_per_molecule_nt for s in samples}
    # spike-in: constant molecules in every fraction
    rows["SIRV1"] = {s: spikein_molecules * spikein_length * depth[s] * reads_per_molecule_nt for s in samples}

    cols = {}
    for s in samples:
        for r in replicates:
            cols[f"{s}_{r}"] = [rows[g][s] for g in rows]
    counts = pd.DataFrame(cols, index=list(rows.keys()))
    if noise:
        counts = pd.DataFrame(rng.poisson(counts.values), index=counts.index, columns=counts.columns)
    counts = counts.reset_index(names=["gene_ID"])

    len_df = pd.DataFrame(
        {
            "gene_ID": gene_ids + ["SIRV1"],
            "transcript_ID": [f"ENST{i:011d}" for i in range(n)] + ["SIRV1"],
            "transcript_length": list(lengths) + [spikein_length],
            "protein_coding": [True] * n + [False],
        }
    )
    return counts, len_df
