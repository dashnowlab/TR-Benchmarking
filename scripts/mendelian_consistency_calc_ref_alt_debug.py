#!/usr/bin/env python3
import sys
import pathlib
import typing as ty
from itertools import product

import plotly.express as px
import cyvcf2


# -----------------------------
# Utilities: chromosome ordering
# -----------------------------
def chrom_key(chrom: str) -> ty.Tuple[int, int, str]:
    """
    Sort key for CHROM that roughly matches genome order:
    1..22, X, Y, M/MT, then others lexicographically.
    Handles chr prefix.
    """
    c = chrom
    if c.startswith("chr"):
        c = c[3:]
    c_upper = c.upper()

    if c_upper.isdigit():
        return (0, int(c_upper), "")
    if c_upper == "X":
        return (1, 23, "")
    if c_upper == "Y":
        return (1, 24, "")
    if c_upper in ("M", "MT"):
        return (1, 25, "")
    return (2, 10**9, c_upper)


def locus_key(
    v: cyvcf2.Variant,
    *,
    merge_key: str = "pos",
) -> ty.Union[ty.Tuple[ty.Tuple[int, int, str], int], str]:
    """
    merge_key="pos" -> (chrom_order, POS)   [original behavior]
    merge_key="id"  -> ID only             [NEW behavior]
    """
    if merge_key == "id":
        return v.ID if v.ID is not None and v.ID != "" else "."
    return (chrom_key(v.CHROM), int(v.POS))


# -----------------------------
# Distance functions (NO Levenshtein)
# -----------------------------
def minkowski(a: ty.List[float], b: ty.List[float], power: int = 1) -> float:
    """power=1 -> Manhattan, power=2 -> Euclidean"""
    return (sum(abs(a[i] - b[i]) ** power for i in range(len(a))) ** (1.0 / power))


def minkowski_sum_only(a: ty.List[float], b: ty.List[float], power: int = 1) -> float:
    """Return the sum(abs(diff)^power) without the final root."""
    return sum(abs(a[i] - b[i]) ** power for i in range(len(a)))


# -----------------------------
# VCF parsing helpers (REF/ALT/GT -> allele sequences)
# -----------------------------
def fmt_gt(v: cyvcf2.Variant) -> str:
    """Return GT like '0/1' (unphased) or '0|1' (phased) if present; else '.'."""
    try:
        g = v.genotypes[0]  # [a, b, phased_bool]
        a, b = g[0], g[1]
        phased = bool(g[2]) if len(g) >= 3 else False
        sep = "|" if phased else "/"
        if a is None or b is None or a < 0 or b < 0:
            return "."
        return f"{a}{sep}{b}"
    except Exception:
        return "."


def allele_seq_from_index(ref: str, alts: ty.List[str], idx: int) -> ty.Optional[str]:
    if idx == 0:
        return ref
    if idx < 0:
        return None
    j = idx - 1
    if j < 0 or j >= len(alts):
        return None
    alt = alts[j]
    # If symbolic ALT like <STR>, can't compute sequence distances
    if alt.startswith("<") and alt.endswith(">"):
        return None
    return alt


def get_allele_seqs(v: cyvcf2.Variant) -> ty.Tuple[ty.Optional[str], ty.Optional[str]]:
    """
    Use REF/ALT + GT to reconstruct the two allele sequences for the first sample.
    Returns (a1, a2) where each can be None if unavailable/symbolic/missing.
    """
    ref = v.REF
    alts = list(v.ALT) if v.ALT is not None else []
    try:
        g = v.genotypes[0]
        a_idx, b_idx = int(g[0]), int(g[1])
    except Exception:
        return (None, None)

    a1 = allele_seq_from_index(ref, alts, a_idx)
    a2 = allele_seq_from_index(ref, alts, b_idx)
    return (a1, a2)


def length_features(ref: str, a1: str, a2: str) -> ty.List[float]:
    """
    Sequence-derived "AL/AP-like" features in bp-space:
      [L1, L2, dL1, dL2] where dL = L - Lref
    """
    Lref = len(ref)
    L1, L2 = len(a1), len(a2)
    return [float(L1), float(L2), float(L1 - Lref), float(L2 - Lref)]


# -----------------------------
# Trio inheritance minimization ("dropdown" selection)
# -----------------------------
def parental_combos(
    mom_a: ty.Tuple[str, str],
    dad_a: ty.Tuple[str, str],
) -> ty.List[ty.Tuple[str, str, str]]:
    """
    4 transmissions: choose one allele from mom and one from dad.
    Returns list of (mom_allele, dad_allele, label) where label is m{0|1}d{0|1}.
    """
    out = []
    for mi, di in product([0, 1], [0, 1]):
        out.append((mom_a[mi], dad_a[di], f"m{mi}d{di}"))
    return out


def min_inheritance_len_distance(
    mom: ty.Tuple[str, str],
    dad: ty.Tuple[str, str],
    kid: ty.Tuple[str, str],
    ref: str,
    power: int = 1,
) -> ty.Tuple[float, str, ty.Tuple[str, str], str]:
    """
    Distance-only selection of best parental transmission and best kid ordering.

    Returns:
      (len_dist, parent_ht_label, kid_ht_pair, kid_ht_ordered)
    """
    # Try both kid orderings
    kid_orders = [
        (kid[0], kid[1], "as_is"),
        (kid[1], kid[0], "flipped"),
    ]

    best = None  # (dist, parent_label, kid_pair, order_tag)

    for k1, k2, order_tag in kid_orders:
        kid_vec = length_features(ref, k1, k2)
        for m_a, d_a, lab in parental_combos(mom, dad):
            p_vec = length_features(ref, m_a, d_a)
            dist = minkowski(p_vec, kid_vec, power=power)
            cand = (dist, lab, (k1, k2), order_tag)
            if best is None or cand[0] < best[0]:
                best = cand

    assert best is not None
    best_dist, best_lab, best_kid_pair, _best_order_tag = best

    # Ordering determinability (same idea as your old code):
    # if best achievable distance differs between the two kid orderings -> ordered=T
    d_as_is = None
    kid_vec_as = length_features(ref, kid[0], kid[1])
    for m_a, d_a, _lab in parental_combos(mom, dad):
        p_vec = length_features(ref, m_a, d_a)
        dist = minkowski(p_vec, kid_vec_as, power=power)
        d_as_is = dist if d_as_is is None else min(d_as_is, dist)

    d_flip = None
    kid_vec_fl = length_features(ref, kid[1], kid[0])
    for m_a, d_a, _lab in parental_combos(mom, dad):
        p_vec = length_features(ref, m_a, d_a)
        dist = minkowski(p_vec, kid_vec_fl, power=power)
        d_flip = dist if d_flip is not None else dist if d_flip is None else min(d_flip, dist)

    kid_ht_ordered = "T" if d_as_is != d_flip else "F"

    return best_dist, best_lab, best_kid_pair, kid_ht_ordered


# -----------------------------
# Merge by key across VCFs
# -----------------------------
def advance(vit):
    try:
        return next(vit)
    except StopIteration:
        return None


def multiway_merge_by_pos(
    vcfs: ty.List[cyvcf2.VCF],
    *,
    merge_key: str = "pos",
) -> ty.Iterator[ty.Tuple[ty.Optional[cyvcf2.Variant], ...]]:
    """
    Stream-merge variants across multiple VCFs by:
      - merge_key="pos": (CHROM, POS) [original]
      - merge_key="id" : ID only
    Assumes each VCF is sorted by the chosen key.
    """
    iters = [iter(v) for v in vcfs]
    cur = [advance(it) for it in iters]

    while any(v is not None for v in cur):
        keys = [locus_key(v, merge_key=merge_key) if v is not None else None for v in cur]
        min_key = min(k for k in keys if k is not None)

        out = []
        for i, v in enumerate(cur):
            if v is not None and locus_key(v, merge_key=merge_key) == min_key:
                out.append(v)
                cur[i] = advance(iters[i])
            else:
                out.append(None)
        yield tuple(out)



def _alts_str(v: cyvcf2.Variant) -> str:
    alts = list(v.ALT) if v.ALT is not None else []
    return ",".join(alts) if alts else "."


def _pick_transmitted_alleles_from_label(
    mom_a: ty.Tuple[str, str],
    dad_a: ty.Tuple[str, str],
    parent_ht: str,
) -> ty.Tuple[str, str]:
    # parent_ht like "m0d1"
    mi = int(parent_ht[1])
    di = int(parent_ht[3])
    return mom_a[mi], dad_a[di]


# -----------------------------
# Main logic
# -----------------------------
def run(
    mom_vcf: pathlib.Path,
    dad_vcf: pathlib.Path,
    kid_vcfs: ty.List[pathlib.Path],
    *,
    power: int = 1,
    output_prefix: str = "seq-mc-",
    exclude_chroms: ty.List[str] = None,
    merge_key: str = "pos",
):
    """
    Sequence-length-feature Mendelian check using only REF/ALT/GT (NO INFO/FORMAT).
    Consistency rule: len_dist == 0  -> consistent, else inconsistent.
    """
    if exclude_chroms is None:
        exclude_chroms = ["chrX", "chrY"]

    vcfs = [cyvcf2.VCF(str(mom_vcf)), cyvcf2.VCF(str(dad_vcf))]
    vcfs.extend([cyvcf2.VCF(str(p)) for p in kid_vcfs])

    kid_ids = [
        vcf.samples[0] if vcf.samples else f"kid{i}"
        for i, vcf in enumerate(vcfs[2:])
    ]

    out_txt = output_prefix + "dists.txt"
    out_html = output_prefix + "dists.html"

    fh = open(out_txt, "w")

    header = (
        "#chrom\tpos\tkid_id\tkid_ID\t"
        "mom_GT\tdad_GT\tkid_GT\t"
        "mom_REF\tmom_ALT\t"
        "dad_REF\tdad_ALT\t"
        "kid_REF\tkid_ALT\t"
        "mom_a1_seq\tmom_a2_seq\t"
        "dad_a1_seq\tdad_a2_seq\t"
        "kid_a1_seq\tkid_a2_seq\t"
        "Lref\t"
        "mom_len\tmom_dlen\t"
        "dad_len\tdad_dlen\t"
        "kid_len\tkid_dlen\t"
        "parent_ht\tmom_tx_seq\tdad_tx_seq\t"
        "kid_ht\tkid_ht_ordered\t"
        "absdiff_L1\tabsdiff_L2\tabsdiff_dL1\tabsdiff_dL2\t"
        "distance_power_sum\tlen_dist\tmendelian"
    )
    print(header, file=fh)

    # For summary / histogram
    len_dists_all: ty.List[float] = []
    total_scored = 0
    total_consistent = 0

    for recs in multiway_merge_by_pos(vcfs, merge_key=merge_key):
        mom, dad = recs[0], recs[1]
        kids = recs[2:]

        if mom is None or dad is None:
            continue
        if mom.CHROM in exclude_chroms:
            continue

        chrom = mom.CHROM
        pos = mom.POS
        ref = mom.REF  # anchor REF from mom
        Lref = len(ref)

        mom_a = get_allele_seqs(mom)
        dad_a = get_allele_seqs(dad)

        # Require mom & dad alleles usable
        if (
            mom_a[0] is None or mom_a[1] is None
            or dad_a[0] is None or dad_a[1] is None
        ):
            continue

        mom_vec = length_features(ref, mom_a[0], mom_a[1])
        dad_vec = length_features(ref, dad_a[0], dad_a[1])
        mom_len = f"{int(mom_vec[0])},{int(mom_vec[1])}"
        mom_dlen = f"{int(mom_vec[2])},{int(mom_vec[3])}"
        dad_len = f"{int(dad_vec[0])},{int(dad_vec[1])}"
        dad_dlen = f"{int(dad_vec[2])},{int(dad_vec[3])}"

        for i, kid in enumerate(kids):
            if kid is None:
                continue

            kid_a = get_allele_seqs(kid)
            if kid_a[0] is None or kid_a[1] is None:
                continue

            kid_vec = length_features(ref, kid_a[0], kid_a[1])
            kid_len = f"{int(kid_vec[0])},{int(kid_vec[1])}"
            kid_dlen = f"{int(kid_vec[2])},{int(kid_vec[3])}"

            len_dist, parent_ht, kid_ht_pair, kid_ht_ordered = min_inheritance_len_distance(
                (mom_a[0], mom_a[1]),
                (dad_a[0], dad_a[1]),
                (kid_a[0], kid_a[1]),
                ref=ref,
                power=power,
            )

            # Reconstruct the chosen transmitted alleles based on parent_ht
            mom_tx, dad_tx = _pick_transmitted_alleles_from_label(
                (mom_a[0], mom_a[1]),
                (dad_a[0], dad_a[1]),
                parent_ht,
            )

            # Parent vector and chosen kid vector for the reported best pairing
            p_vec = length_features(ref, mom_tx, dad_tx)
            k_vec = length_features(ref, kid_ht_pair[0], kid_ht_pair[1])

            absdiff = [abs(p_vec[j] - k_vec[j]) for j in range(4)]
            dist_sum = minkowski_sum_only(p_vec, k_vec, power=power)

            mendelian = "T" if len_dist == 0 else "F"

            total_scored += 1
            if mendelian == "T":
                total_consistent += 1
            len_dists_all.append(len_dist)

            # kid_ht: keep compact as allele lengths of the chosen ordering
            kid_ht = f"{len(kid_ht_pair[0])},{len(kid_ht_pair[1])}"

            # kid_ID from VCF ID column (no logic use)
            kid_ID = kid.ID if kid.ID is not None and kid.ID != "" else "."

            line = (
                f"{chrom}\t{pos}\t{kid_ids[i]}\t{kid_ID}\t"
                f"{fmt_gt(mom)}\t{fmt_gt(dad)}\t{fmt_gt(kid)}\t"
                f"{mom.REF}\t{_alts_str(mom)}\t"
                f"{dad.REF}\t{_alts_str(dad)}\t"
                f"{kid.REF}\t{_alts_str(kid)}\t"
                f"{mom_a[0]}\t{mom_a[1]}\t"
                f"{dad_a[0]}\t{dad_a[1]}\t"
                f"{kid_a[0]}\t{kid_a[1]}\t"
                f"{Lref}\t"
                f"{mom_len}\t{mom_dlen}\t"
                f"{dad_len}\t{dad_dlen}\t"
                f"{kid_len}\t{kid_dlen}\t"
                f"{parent_ht}\t{mom_tx}\t{dad_tx}\t"
                f"{kid_ht}\t{kid_ht_ordered}\t"
                f"{int(absdiff[0])}\t{int(absdiff[1])}\t{int(absdiff[2])}\t{int(absdiff[3])}\t"
                f"{dist_sum}\t{len_dist}\t{mendelian}"
            )
            print(line, file=fh)

    fh.close()

    # Summary like original (percentage)
    if total_scored > 0:
        pct = (total_consistent / total_scored) * 100.0
    else:
        pct = 0.0
    print(
        f"[summary] scored_loci={total_scored} consistent={total_consistent} "
        f"pct={pct:.2f}%",
        file=sys.stderr,
    )

    # Histogram
    if len(len_dists_all) > 0:
        mname = "manhattan" if power == 1 else ("euclidean" if power == 2 else f"minkowski(p={power})")
        fig = px.histogram(x=len_dists_all)
        fig.update_layout(
            xaxis_title=f"{mname} distance using length features [L1,L2,dL1,dL2] (min transmission)",
            yaxis_title="count",
        )
        fig.write_html(out_html)
        print(f"wrote {out_html}", file=sys.stderr)

    print(f"wrote {out_txt}", file=sys.stderr)


def parse_args(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Sequence-length Mendelian consistency using REF/ALT + GT only (no INFO/FORMAT dependence)."
    )
    p.add_argument("--mom", required=True, type=pathlib.Path, help="Mom VCF (.vcf/.vcf.gz)")
    p.add_argument("--dad", required=True, type=pathlib.Path, help="Dad VCF (.vcf/.vcf.gz)")
    p.add_argument(
        "--kids",
        required=True,
        action="append",
        type=pathlib.Path,
        help="Kid VCF (.vcf/.vcf.gz). Repeat --kids for multiple kids.",
    )
    p.add_argument(
        "--power",
        type=int,
        default=1,
        help="1=Manhattan, 2=Euclidean for length-feature distance (default: 1).",
    )
    p.add_argument(
        "--out-prefix",
        default="seq-mc-",
        help="Output prefix (default: seq-mc-).",
    )
    p.add_argument(
        "--exclude-chroms",
        nargs="*",
        default=["chrX", "chrY"],
        help="Chromosomes to skip (default: chrX chrY).",
    )
    p.add_argument(
        "--merge-key",
        choices=["pos", "id"],
        default="pos",
        help="How to merge loci across VCFs: 'pos' (CHROM+POS, default) or 'id' (CHROM+ID).",
    )
    return p.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    run(
        args.mom,
        args.dad,
        args.kids,
        power=args.power,
        output_prefix=args.out_prefix,
        exclude_chroms=args.exclude_chroms,
        merge_key=args.merge_key,
    )
