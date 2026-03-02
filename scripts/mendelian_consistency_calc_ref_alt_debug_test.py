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
# Distance functions
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


def length_features(a1: str, a2: str) -> ty.List[float]:
    """
    FIX 1: Removed redundant dL features.
    The original [L1, L2, dL1, dL2] was redundant because dL = L - Lref
    (constant per locus), so dL added zero discriminative information and
    inflated distances by a factor of 2.
    Now returns just [L1, L2] in bp-space, analogous to the AL tag in TRGT.
    """
    return [float(len(a1)), float(len(a2))]


# -----------------------------
# Trio inheritance minimization
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
    power: int = 1,
    consistency_threshold: float = 0.0,
) -> ty.Tuple[float, str, ty.Tuple[str, str], str]:
    """
    Distance-based selection of best parental transmission and best kid ordering.

    FIX 2: Removed ref parameter — it was only used by the now-removed dL features.

    FIX 3: Fixed d_flip initialization bug.
    The original code had:
        d_flip = dist if d_flip is not None else dist if d_flip is None else min(d_flip, dist)
    The min() branch was unreachable (both conditions cover all cases), so d_flip
    always ended up as the LAST distance seen, not the MINIMUM. Fixed to:
        d_flip = dist if d_flip is None else min(d_flip, dist)

    Returns:
      (len_dist, parent_ht_label, kid_ht_pair, kid_ht_ordered)
    """
    # Try both kid orderings and all parental combinations
    kid_orders = [
        (kid[0], kid[1], "as_is"),
        (kid[1], kid[0], "flipped"),
    ]

    best = None  # (dist, parent_label, kid_pair, order_tag)

    for k1, k2, order_tag in kid_orders:
        kid_vec = length_features(k1, k2)
        for m_a, d_a, lab in parental_combos(mom, dad):
            p_vec = length_features(m_a, d_a)
            dist = minkowski(p_vec, kid_vec, power=power)
            cand = (dist, lab, (k1, k2), order_tag)
            if best is None or cand[0] < best[0]:
                best = cand

    assert best is not None
    best_dist, best_lab, best_kid_pair, _best_order_tag = best

    # Compute minimum achievable distance for each kid ordering independently
    # to determine whether allele phase can be resolved (kid_ht_ordered).
    d_as_is = None
    kid_vec_as = length_features(kid[0], kid[1])
    for m_a, d_a, _lab in parental_combos(mom, dad):
        p_vec = length_features(m_a, d_a)
        dist = minkowski(p_vec, kid_vec_as, power=power)
        # FIX 3 (as_is side — was already correct in original, kept for symmetry):
        d_as_is = dist if d_as_is is None else min(d_as_is, dist)

    # FIX 3: d_flip bug fixed here
    d_flip = None
    kid_vec_fl = length_features(kid[1], kid[0])
    for m_a, d_a, _lab in parental_combos(mom, dad):
        p_vec = length_features(m_a, d_a)
        dist = minkowski(p_vec, kid_vec_fl, power=power)
        d_flip = dist if d_flip is None else min(d_flip, dist)  # FIXED

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
    consistency_threshold: float = 0.0,
):
    """
    Sequence-length-feature Mendelian check using only REF/ALT/GT (NO INFO/FORMAT).

    FIX 4: consistency_threshold is now a parameter (default 0.0 = exact match).
    For long-read tools (Medaka, STRdust, etc.) that may produce slightly varying
    allele lengths for the same true allele, consider passing a small threshold
    (e.g. --consistency-threshold 2) to avoid overcounting inconsistency.

    FIX 5: symbolic ALT skips are now counted and reported per-kid so the
    denominator difference across tools is visible and not silently hidden.
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
        "mom_L1\tmom_L2\t"
        "dad_L1\tdad_L2\t"
        "kid_L1\tkid_L2\t"
        "parent_ht\tmom_tx_seq\tdad_tx_seq\t"
        "kid_ht\tkid_ht_ordered\t"
        "absdiff_L1\tabsdiff_L2\t"
        "distance_power_sum\tlen_dist\tmendelian"
    )
    print(header, file=fh)

    # Per-kid tracking: FIX 5
    n_kids = len(kid_vcfs)
    total_scored = [0] * n_kids
    total_consistent = [0] * n_kids
    total_skipped_symbolic = [0] * n_kids   # FIX 5: track symbolic ALT skips
    total_skipped_missing = [0] * n_kids    # FIX 5: track missing GT skips

    len_dists_all: ty.List[float] = []

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

        mom_a = get_allele_seqs(mom)
        dad_a = get_allele_seqs(dad)

        # Require mom & dad alleles usable
        if (
            mom_a[0] is None or mom_a[1] is None
            or dad_a[0] is None or dad_a[1] is None
        ):
            continue

        mom_vec = length_features(mom_a[0], mom_a[1])
        dad_vec = length_features(dad_a[0], dad_a[1])

        for i, kid in enumerate(kids):
            if kid is None:
                continue

            kid_a = get_allele_seqs(kid)

            # FIX 5: distinguish why a locus is skipped
            if kid_a[0] is None or kid_a[1] is None:
                # Check if it's symbolic vs missing GT
                try:
                    g = kid.genotypes[0]
                    a_idx, b_idx = int(g[0]), int(g[1])
                    alts = list(kid.ALT) if kid.ALT is not None else []
                    # If indices are valid but returned None -> symbolic ALT
                    if a_idx >= 0 and b_idx >= 0:
                        total_skipped_symbolic[i] += 1
                    else:
                        total_skipped_missing[i] += 1
                except Exception:
                    total_skipped_missing[i] += 1
                continue

            kid_vec = length_features(kid_a[0], kid_a[1])

            # FIX 2: removed ref= argument (no longer needed after FIX 1)
            len_dist, parent_ht, kid_ht_pair, kid_ht_ordered = min_inheritance_len_distance(
                (mom_a[0], mom_a[1]),
                (dad_a[0], dad_a[1]),
                (kid_a[0], kid_a[1]),
                power=power,
                consistency_threshold=consistency_threshold,
            )

            # Reconstruct the chosen transmitted alleles based on parent_ht
            mom_tx, dad_tx = _pick_transmitted_alleles_from_label(
                (mom_a[0], mom_a[1]),
                (dad_a[0], dad_a[1]),
                parent_ht,
            )

            # Parent vector and chosen kid vector for the reported best pairing
            p_vec = length_features(mom_tx, dad_tx)
            k_vec = length_features(kid_ht_pair[0], kid_ht_pair[1])

            absdiff = [abs(p_vec[j] - k_vec[j]) for j in range(2)]
            dist_sum = minkowski_sum_only(p_vec, k_vec, power=power)

            # FIX 4: use configurable threshold instead of hard-coded 0
            mendelian = "T" if len_dist <= consistency_threshold else "F"

            total_scored[i] += 1
            if mendelian == "T":
                total_consistent[i] += 1
            len_dists_all.append(len_dist)

            # kid_ht: allele lengths of the chosen ordering
            kid_ht = f"{len(kid_ht_pair[0])},{len(kid_ht_pair[1])}"

            # kid_ID from VCF ID column (no logic use)
            kid_ID = kid.ID if kid.ID is not None and kid.ID != "" else "."

            # FIX 6: corrected kid.REF == kid.REF typo -> mom.REF == kid.REF
            if mom.REF != kid.REF:
                sys.stderr.write(
                    f"Warning: REF mismatch at {chrom}:{pos} between mom and kid {kid_ids[i]}. "
                    f"mom={mom.REF} kid={kid.REF}\n"
                )

            line = (
                f"{chrom}\t{pos}\t{kid_ids[i]}\t{kid_ID}\t"
                f"{fmt_gt(mom)}\t{fmt_gt(dad)}\t{fmt_gt(kid)}\t"
                f"{mom.REF}\t{_alts_str(mom)}\t"
                f"{dad.REF}\t{_alts_str(dad)}\t"
                f"{kid.REF}\t{_alts_str(kid)}\t"
                f"{mom_a[0]}\t{mom_a[1]}\t"
                f"{dad_a[0]}\t{dad_a[1]}\t"
                f"{kid_a[0]}\t{kid_a[1]}\t"
                f"{int(mom_vec[0])}\t{int(mom_vec[1])}\t"
                f"{int(dad_vec[0])}\t{int(dad_vec[1])}\t"
                f"{int(kid_vec[0])}\t{int(kid_vec[1])}\t"
                f"{parent_ht}\t{mom_tx}\t{dad_tx}\t"
                f"{kid_ht}\t{kid_ht_ordered}\t"
                f"{int(absdiff[0])}\t{int(absdiff[1])}\t"
                f"{dist_sum}\t{len_dist}\t{mendelian}"
            )
            print(line, file=fh)

    fh.close()

    # Per-kid summary (FIX 5: includes skipped counts)
    print("[summary]", file=sys.stderr)
    for i, kid_id in enumerate(kid_ids):
        scored = total_scored[i]
        consistent = total_consistent[i]
        skipped_sym = total_skipped_symbolic[i]
        skipped_mis = total_skipped_missing[i]
        pct = (consistent / scored * 100.0) if scored > 0 else 0.0
        print(
            f"  kid={kid_id} scored={scored} consistent={consistent} "
            f"pct={pct:.2f}% skipped_symbolic={skipped_sym} skipped_missing={skipped_mis}",
            file=sys.stderr,
        )

    # Histogram
    if len(len_dists_all) > 0:
        mname = "manhattan" if power == 1 else ("euclidean" if power == 2 else f"minkowski(p={power})")
        fig = px.histogram(x=len_dists_all)
        fig.update_layout(
            xaxis_title=f"{mname} distance using length features [L1,L2] (min transmission)",
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
    # FIX 4: expose threshold as CLI argument
    p.add_argument(
        "--consistency-threshold",
        type=float,
        default=0.0,
        help=(
            "Maximum allele-length distance to call a locus Mendelian-consistent "
            "(default: 0.0 = exact match). For long-read tools with noisy length "
            "estimates consider a small value like 2."
        ),
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
        consistency_threshold=args.consistency_threshold,
    )