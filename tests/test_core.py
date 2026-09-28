"""Unit tests on small synthetic inputs (no downloads needed). Run with ``pytest``."""

import itertools
import random

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from phosphokin.align import BLOSUM62, _encode, align_position_map
from phosphokin.enrichment import bh_adjust, classify_sites, kinase_enrichment
from phosphokin.kinase_atlas import POSITIONS, RESIDUES, KinomeAtlas
from phosphokin.orthology import map_sites_to_human, orthologs_by_symbol
from phosphokin.pathways import overrepresentation
from phosphokin.peptides import parse_probability_sequence, peptides_to_sites
from phosphokin.sequences import center_window, site_window


# ----------------------------------------------------------------- sequences
def test_site_window_pads_termini():
    seq = "MSTPKRAAG"
    assert site_window(seq, 2) == "______MSTPKRAAG"
    assert site_window(seq, 9, flank=3) == "RAAG___"
    assert center_window("A" * 15 + "S" + "B" * 15) == "AAAAAAASBBBBBBB"


# ------------------------------------------------------------------ peptides
def test_parse_probability_sequence():
    pep, sites = parse_probability_sequence("AAS(0.98)PT(0.02)KY(0)R")
    assert pep == "AASPTKYR"
    assert sites == [("S", 3, 0.98), ("T", 5, 0.02)]


def test_peptides_to_sites_maps_and_filters():
    proteins = pd.DataFrame({"accession": ["P1"], "gene": ["Abc"], "reviewed": [True],
                             "sequence": ["MKKLLAASPTKRRE"]})
    peps = pd.DataFrame({"accession": ["P1", "P1"],
                         "annotated_peptide": ["AAS(0.98)PT(0.02)K", "QQQS(1)K"]})
    out = peptides_to_sites(peps, proteins)
    assert out["site"].tolist() == ["S8"]
    assert out.attrs["unmapped"] == [("P1", "QQQSK")]


# ----------------------------------------------------------------- alignment
def _reference_score(a, b, go, ge):
    A, B = _encode(a), _encode(b)
    n, m, neg = len(A), len(B), -10 ** 9
    H = [[neg] * (m + 1) for _ in range(n + 1)]
    X = [[neg] * (m + 1) for _ in range(n + 1)]
    Y = [[neg] * (m + 1) for _ in range(n + 1)]
    H[0][0] = 0
    for j in range(1, m + 1):
        H[0][j] = Y[0][j] = -(go + (j - 1) * ge)
    for i in range(1, n + 1):
        H[i][0] = X[i][0] = -(go + (i - 1) * ge)
    for i, j in itertools.product(range(1, n + 1), range(1, m + 1)):
        X[i][j] = max(H[i - 1][j] - go, X[i - 1][j] - ge)
        Y[i][j] = max(H[i][j - 1] - go, Y[i][j - 1] - ge)
        H[i][j] = max(H[i - 1][j - 1] + int(BLOSUM62[A[i - 1], B[j - 1]]), X[i][j], Y[i][j])
    return H[n][m]


def test_blosum62_is_symmetric():
    assert (BLOSUM62 == BLOSUM62.T).all()
    assert BLOSUM62[_encode("W")[0], _encode("W")[0]] == 11


def test_alignment_matches_brute_force():
    rng = random.Random(0)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    for _ in range(100):
        a = "".join(rng.choice(aa) for _ in range(rng.randint(1, 20)))
        b = "".join(rng.choice(aa) for _ in range(rng.randint(1, 20)))
        _, score = align_position_map(a, b)
        assert score == _reference_score(a, b, 11, 1)


def test_alignment_maps_across_an_insertion():
    a = "MSTAAKRPLSEEGHWY"
    b = "MSTAAKGGGGRPLSEEGHWY"          # 4-residue insertion in b
    pos_map, _ = align_position_map(a, b, gap_open=5, gap_extend=1)
    assert pos_map[10] == 14            # S10 in a -> S14 in b
    assert b[pos_map[10] - 1] == "S"


# ----------------------------------------------------------------- orthology
def test_map_sites_to_human_alignment_and_window_agree():
    mouse_seq = "MSDNQKRPLSPEGHAVKTEEDLLSGSRPASPT"
    human_seq = "MSENQKRAPLSPEGHAVKTEEDLLSGSRPASPT"   # 1 substitution + 1 insertion
    mouse = pd.DataFrame({"accession": ["Q1"], "gene": ["Abc1"], "reviewed": [True], "sequence": [mouse_seq]})
    human = pd.DataFrame({"accession": ["H1"], "gene": ["ABC1"], "reviewed": [True], "sequence": [human_seq]})
    sites = pd.DataFrame({"accession": ["Q1", "Q1"], "gene": ["Abc1", "Abc1"],
                          "residue": ["S", "S"], "position": [10, 24]})
    orth = orthologs_by_symbol(sites["gene"], human)
    a = map_sites_to_human(sites, mouse, human, orth, method="alignment")
    w = map_sites_to_human(sites, mouse, human, orth, method="window")
    assert a["human_site"].tolist() == ["S11", "S25"]
    assert w["human_site"].tolist() == a["human_site"].tolist()


def test_map_sites_relocates_shifted_positions():
    seq = "MKKLLAASPTKRRESPEQLVVGT"
    src = pd.DataFrame({"accession": ["Q1"], "gene": ["Abc1"], "reviewed": [True], "sequence": [seq]})
    hum = pd.DataFrame({"accession": ["H1"], "gene": ["ABC1"], "reviewed": [True], "sequence": [seq]})
    # position is off by 3 (different sequence version), but the measured window is correct
    sites = pd.DataFrame({"accession": ["Q1"], "gene": ["Abc1"], "residue": ["S"], "position": [5],
                          "window": [site_window(seq, 8)]})
    out = map_sites_to_human(sites, src, hum, orthologs_by_symbol(["Abc1"], hum))
    assert int(out["source_position"].iat[0]) == 8 and out["human_site"].iat[0] == "S8"


# ------------------------------------------------------------- kinase atlas
def _toy_atlas():
    kinases = ["BASO", "ACIDO"]          # likes R at -3 / likes E at +1
    cols = [f"{p}{aa}" for p in POSITIONS for aa in RESIDUES]
    m = pd.DataFrame(1.0, index=kinases, columns=cols)
    m.loc["BASO", "-3R"] = 8.0
    m.loc["ACIDO", "1E"] = 8.0
    return KinomeAtlas.from_tables(m)


def test_scores_follow_motifs():
    atlas = _toy_atlas()
    s = atlas.score(["AAAARAASAAAAAAA", "AAAAAAASEAAAAAA", "_______SAAAAAAA"])
    assert np.allclose(s, [[3, 0], [0, 3], [0, 0]])   # log2(8) = 3; padding is neutral


def test_percentiles_and_ranks_from_reference():
    atlas = _toy_atlas()
    rng = np.random.default_rng(1)
    windows = ["".join(rng.choice(list("ARDEGLKS"), 15)[:7]) + "S" + "".join(rng.choice(list("ARDEGLKS"), 7))
               for _ in range(400)]
    scores = atlas.score(windows)
    ref = pd.DataFrame({"SITE_+/-7_AA": windows})
    for k, name in enumerate(atlas.kinases):     # published percentiles = empirical CDF of scores
        ref[f"{name}_percentile"] = pd.Series(scores[:, k]).rank(pct=True).to_numpy() * 100
    atlas.set_reference(ref)
    pct = atlas.percentile(windows)
    assert np.allclose(pct, ref[["BASO_percentile", "ACIDO_percentile"]].to_numpy())
    ranks = atlas.rank(atlas.percentile(["AAAARAASAAAAAAA"]))
    assert ranks.tolist() == [[1, 2]]


# ---------------------------------------------------------------- enrichment
def test_bh_adjust_known_values():
    p = [0.01, 0.04, 0.03, 0.005]
    assert np.allclose(bh_adjust(p), [0.02, 0.04, 0.04, 0.02])
    assert np.isnan(bh_adjust([np.nan, 0.5])[0])


def test_classify_sites():
    lab = classify_sites([2.5, -3, 0.2, 2.5], pvalue=[0.01, 0.01, 0.01, 0.2], fc_threshold=1)
    assert lab.tolist() == ["up", "down", "unregulated", "unregulated"]


def test_kinase_enrichment_matches_fisher():
    status = pd.DataFrame({"site_id": [f"s{i}" for i in range(40)],
                           "status": ["up"] * 10 + ["down"] * 10 + ["unregulated"] * 20})
    # kinase K favours 8/10 up sites and 2/20 unregulated sites; every site has one favoured kinase
    fav = pd.DataFrame({"site_id": [f"s{i}" for i in range(40)],
                        "kinase": ["K" if (i < 8 or 20 <= i < 22) else "Z" for i in range(40)]})
    res = kinase_enrichment(status, fav).set_index("kinase")
    p = fisher_exact([[8, 2], [2, 18]], alternative="greater")[1]
    assert np.isclose(res.at["K", "p_up"], p)
    assert np.isclose(res.at["K", "log2_ff_up"], np.log2((8 / 10) / (2 / 20)))
    assert res.at["K", "direction"] == "up"


def test_overrepresentation():
    sets = {"A": {f"g{i}" for i in range(10)}, "B": {f"g{i}" for i in range(50, 60)}}
    res = overrepresentation([f"g{i}" for i in range(8)], [f"g{i}" for i in range(100)], sets)
    assert res.iloc[0]["gene_set"] == "A" and res.iloc[0]["overlap"] == 8
