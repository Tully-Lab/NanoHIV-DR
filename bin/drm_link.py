#!/usr/bin/env python3
"""
DRM-LINK: read-backed DRM co-occurrence / within-amplicon linkage (prototype).

Reports whether pairs of HIV-1 drug-resistance mutations (DRMs) are carried on
the same sequenced molecule (within the amplicon) or segregate across molecules.

IMPORTANT (see the DRM-LINK spec):
  * A read is a single amplified molecule, NOT a viral genome.
  * Because the input is a PCR amplicon, co-occurrence on a read reflects
    within-amplicon linkage and may include artefactual PCR recombinants.
    Supply --chimera-floor (measured from an orthogonal-clone control mix) so
    linkage is tested against an independence+chimera-adjusted null, not against
    independence alone.

Codon-aware: amino acids are reconstructed from all three positions of each
codon using the read-to-reference alignment; each locus is scored WT / DRM /
uncertain. Uncertain calls are excluded pairwise (per-pair complete observation).

Usage:
    drm_link.py --bam sample.bam --targets targets.tsv --sample S1 \
        --outdir linkage_out [--chimera-floor 0.005] [--min-mapq 20] \
        [--min-bq 10] [--min-cocov 50] [--min-support 5] [--min-freq 0.01]

    drm_link.py --selftest      # run statistics self-test, no BAM required

targets.tsv (tab-separated, header required):
    gene    position    ref_aa    drm_aa
    RT      65          K         R
    RT      184         M         V
    RT      103         K         N
    RT      181         Y         C
(drm_aa may be several residues, e.g. "IV"; any listed residue counts as DRM.)
"""
import argparse
import csv
import math
import os
import sys
from itertools import combinations

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

# 0-based nucleotide offset of the first base of codon 1 for each gene, in the
# HXB2 pol reference used by NanoHIV-DR (PR: 99 aa, RT: 560 aa, IN: 288 aa).
# Override with --gene-offsets if your reference differs.
GENE_OFFSETS = {"PR": 0, "RT": 297, "IN": 297 + 560 * 3}


def translate(codon):
    codon = (codon or "").upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return None
    return CODON_TABLE.get(codon)


# ----------------------------------------------------------------------
# Pure statistics (unit-testable without a BAM)
# ----------------------------------------------------------------------
def wilson_ci(k, n, z=1.96):
    """Wilson score interval for a binomial proportion."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def odds_ratio(n11, n10, n01, n00):
    """Odds ratio with Haldane-Anscombe (+0.5) correction for zero cells."""
    a, b, c, d = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5
    return (a * d) / (b * c)


def dprime_r2(n11, n10, n01, n00):
    """LD measures treating each locus as biallelic (WT vs DRM).

    n11 = DRM at both A and B; n10 = DRM at A, WT at B; etc.
    Returns (D, D_prime, r2, pA, pB) where pA/pB are DRM frequencies.
    """
    n = n11 + n10 + n01 + n00
    if n == 0:
        return (0.0, 0.0, 0.0, 0.0, 0.0)
    pA = (n11 + n10) / n          # DRM freq at locus A
    pB = (n11 + n01) / n          # DRM freq at locus B
    D = n11 / n - pA * pB
    if D > 0:
        dmax = min(pA * (1 - pB), (1 - pA) * pB)
    else:
        dmax = min(pA * pB, (1 - pA) * (1 - pB))
    dprime = D / dmax if dmax > 0 else 0.0
    denom = pA * (1 - pA) * pB * (1 - pB)
    r2 = (D * D) / denom if denom > 0 else 0.0
    return (D, dprime, r2, pA, pB)


def conditional(n_ab, n_a):
    """P(B | A) = N(A and B) / N(A)."""
    return (n_ab / n_a) if n_a else float("nan")


def excess_over_null(n11, n10, n01, n00, chimera_floor=0.0):
    """Observed double-DRM frequency minus the independence+chimera null.

    The null for genuine linkage is the independence expectation PLUS the
    empirically measured artefactual (chimera) contribution, so that at very
    high depth trivial deviations from independence are not called as linkage.
    """
    n = n11 + n10 + n01 + n00
    if n == 0:
        return 0.0
    obs = n11 / n
    pA = (n11 + n10) / n
    pB = (n11 + n01) / n
    indep = pA * pB
    null = max(indep, indep + chimera_floor)  # floor lifts the null upward
    return obs - null


# ----------------------------------------------------------------------
# Codon-aware BAM parsing (requires pysam)
# ----------------------------------------------------------------------
def read_codon_state(aln, ref_positions, min_bq):
    """Return the amino acid for a read at the 3 given 0-based ref positions,
    or None if any position is unaligned / low quality / indel-disrupted."""
    seq = aln.query_sequence
    quals = aln.query_qualities
    if seq is None:
        return None
    wanted = set(ref_positions)
    got = {}
    for qpos, rpos in aln.get_aligned_pairs(matches_only=False):
        if rpos in wanted:
            if qpos is None:              # deletion in read at this ref pos
                return None
            if quals is not None and quals[qpos] < min_bq:
                return None
            got[rpos] = seq[qpos]
    if len(got) != 3:
        return None
    codon = "".join(got[p] for p in ref_positions)
    return translate(codon)


def parse_bam(bam_path, targets, gene_offsets, min_mapq, min_bq):
    """Return per-read locus states: {read_id: {locus_key: 'WT'|'DRM'|None}}.

    locus_key = f"{gene}{position}{ref_aa}" (e.g. 'RT65K'); the display label
    used downstream is f"{ref_aa}{position}{drm_aa}".
    """
    import pysam  # imported lazily so --selftest needs no pysam
    reads = {}
    bam = pysam.AlignmentFile(bam_path, "rb")
    for t in targets:
        gene, pos = t["gene"], int(t["position"])
        offset = gene_offsets.get(gene)
        if offset is None:
            sys.stderr.write(f"WARNING: no reference offset for gene {gene}; skipping {gene}{pos}\n")
            continue
        codon_start = offset + (pos - 1) * 3
        ref_positions = [codon_start, codon_start + 1, codon_start + 2]
        drm_set = set(t["drm_aa"].upper())
        key = t["_key"]
        for aln in bam.fetch():
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            aa = read_codon_state(aln, ref_positions, min_bq)
            state = None
            if aa is not None:
                state = "DRM" if aa in drm_set else "WT"
            reads.setdefault(aln.query_name, {})[key] = state
        bam.reset()
    bam.close()
    return reads


# ----------------------------------------------------------------------
# Pairwise linkage over parsed reads
# ----------------------------------------------------------------------
def pairwise_table(reads, key_a, key_b):
    """Return (n11, n10, n01, n00) over reads that call BOTH loci (WT/DRM)."""
    n11 = n10 = n01 = n00 = 0
    for states in reads.values():
        a = states.get(key_a)
        b = states.get(key_b)
        if a is None or b is None:
            continue
        if a == "DRM" and b == "DRM":
            n11 += 1
        elif a == "DRM":
            n10 += 1
        elif b == "DRM":
            n01 += 1
        else:
            n00 += 1
    return n11, n10, n01, n00


def assign_tier(n11, cocov, r2, excess, min_support, min_freq):
    if cocov == 0 or n11 == 0:
        return "0_none"
    double_freq = n11 / cocov
    observed = "1_co-occurrence"
    if excess > 0 and r2 >= 0.1 and double_freq >= min_freq:
        observed = "2_linkage"
        if n11 >= min_support:
            observed = "3_high-confidence"
    return observed


def analyse(reads, targets, args):
    label = {t["_key"]: t["_label"] for t in targets}
    rows = []
    for ka, kb in combinations([t["_key"] for t in targets], 2):
        n11, n10, n01, n00 = pairwise_table(reads, ka, kb)
        cocov = n11 + n10 + n01 + n00
        if cocov < args.min_cocov:
            rows.append(dict(a=label[ka], b=label[kb], cocov=cocov,
                             n11=n11, n10=n10, n01=n01, n00=n00,
                             double_freq="", lo="", hi="", OR="", Dprime="",
                             r2="", excess="", pB_given_A="", pA_given_B="",
                             tier="0_insufficient-cocoverage"))
            continue
        D, dprime, r2, pA, pB = dprime_r2(n11, n10, n01, n00)
        lo, hi = wilson_ci(n11, cocov)
        OR = odds_ratio(n11, n10, n01, n00)
        exc = excess_over_null(n11, n10, n01, n00, args.chimera_floor)
        pBgA = conditional(n11, n11 + n10)
        pAgB = conditional(n11, n11 + n01)
        tier = assign_tier(n11, cocov, r2, exc, args.min_support, args.min_freq)
        rows.append(dict(a=label[ka], b=label[kb], cocov=cocov,
                         n11=n11, n10=n10, n01=n01, n00=n00,
                         double_freq=round(n11 / cocov, 4),
                         lo=round(lo, 4), hi=round(hi, 4),
                         OR=round(OR, 3), Dprime=round(dprime, 3), r2=round(r2, 3),
                         excess=round(exc, 4),
                         pB_given_A=round(pBgA, 3) if pBgA == pBgA else "",
                         pA_given_B=round(pAgB, 3) if pAgB == pAgB else "",
                         tier=tier))
    return rows


def write_tsv(rows, path):
    cols = ["a", "b", "cocov", "n11", "n10", "n01", "n00", "double_freq",
            "lo", "hi", "OR", "Dprime", "r2", "excess",
            "pB_given_A", "pA_given_B", "tier"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def load_targets(path):
    targets = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row = {k.strip(): (v or "").strip() for k, v in row.items()}
            row["_key"] = f'{row["gene"]}{row["position"]}{row["ref_aa"]}'
            row["_label"] = f'{row["ref_aa"]}{row["position"]}{row["drm_aa"]}'
            targets.append(row)
    return targets


# ----------------------------------------------------------------------
def selftest():
    # Example 2x2 from the spec (RT65 x RT184): strong positive linkage.
    n11, n10, n01, n00 = 1300, 120, 80, 6500
    D, dprime, r2, pA, pB = dprime_r2(n11, n10, n01, n00)
    OR = odds_ratio(n11, n10, n01, n00)
    lo, hi = wilson_ci(n11, n11 + n10 + n01 + n00)
    print("Self-test (RT65 x RT184 example):")
    print(f"  pA(65R)={pA:.3f}  pB(184V)={pB:.3f}")
    print(f"  double_freq={n11/(n11+n10+n01+n00):.3f}  95% CI=({lo:.3f},{hi:.3f})")
    print(f"  D'={dprime:.3f}  r2={r2:.3f}  OR={OR:.1f}")
    print(f"  P(184V|65R)={conditional(n11, n11+n10):.3f}")
    print(f"  P(65R|184V)={conditional(n11, n11+n01):.3f}")
    print(f"  excess over independence (floor=0)   = {excess_over_null(n11,n10,n01,n00,0.0):.4f}")
    print(f"  excess over independence+floor(0.02) = {excess_over_null(n11,n10,n01,n00,0.02):.4f}")
    # Independence example: no linkage
    D2, dp2, r22, _, _ = dprime_r2(100, 100, 100, 100)
    print(f"  independence control: D'={dp2:.3f} r2={r22:.3f} (expect ~0)")
    ok = r2 > 0.6 and OR > 100 and abs(r22) < 1e-9
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description="DRM-LINK read-backed linkage (prototype)")
    ap.add_argument("--bam")
    ap.add_argument("--targets")
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--outdir", default="linkage_out")
    ap.add_argument("--chimera-floor", type=float, default=0.0,
                    help="Empirical artefactual double-positive frequency (from a control mix)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-bq", type=int, default=10)
    ap.add_argument("--min-cocov", type=int, default=50)
    ap.add_argument("--min-support", type=int, default=5)
    ap.add_argument("--min-freq", type=float, default=0.01)
    ap.add_argument("--gene-offsets",
                    help="Override gene offsets, e.g. 'PR:0,RT:297,IN:1977'")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        return selftest()

    if not args.bam or not args.targets:
        ap.error("--bam and --targets are required (or use --selftest)")

    gene_offsets = dict(GENE_OFFSETS)
    if args.gene_offsets:
        for tok in args.gene_offsets.split(","):
            g, o = tok.split(":")
            gene_offsets[g.strip()] = int(o)

    targets = load_targets(args.targets)
    reads = parse_bam(args.bam, targets, gene_offsets, args.min_mapq, args.min_bq)
    rows = analyse(reads, targets, args)

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"{args.sample}.linkage.tsv")
    write_tsv(rows, out)
    print(f"DRM-LINK: {len(reads)} reads parsed, {len(rows)} locus pairs -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
