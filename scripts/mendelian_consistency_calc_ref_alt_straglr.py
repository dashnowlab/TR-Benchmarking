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
    merge_key="pos" -> (chrom_order, POS)
    merge_key="id"  -> ID only
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
# GT formatting
# -----------------------------
def fmt_gt(v: cyvcf2.Variant) -> str:
    """Return GT like '0/1' or '0|1' if present; else '.'."""
    try:
        g = v.genotypes[0]
        if len(g) < 2:
            return "."
        a, b = g[0], g[1]
        phased = bool(g[2]) if len(g) >= 3 else False
        sep = "|" if phased else "/"
        if a is None or b is None or a < 0 or b < 0:
            return "."
        return f"{a}{sep}{b}"
    except Exception:
        return "."


def _alts_str(v: cyvcf2.Variant) -> str:
    alts = list(v.ALT) if v.ALT is not None else []
    return ",".join(alts) if alts else "."


# -----------------------------
# Straglr: derive allele lengths from INFO (RB, RUL_REF, END) + GT
# -----------------------------
def _parse_info_list(val) -> ty.List[float]:
    """
    cyvcf2 INFO values can be:
      - None
      - int/float
      - string "a,b,c"
      - tuple/list
    Return list[float].
    """
    if val is None:
        return []
    if isinstance(val, (int, float)):
        return [float(val)]
    if isinstance(val, (list, tuple)):
        out = []
        for x in val:
            try:
                out.append(float(x))
            except Exception:
                pass
        return out
    s = str(val).strip()
    if s in ("", "."):
        return []
    out = []
    for tok in s.split(","):
        tok = tok.strip()
        if tok == "":
            continue
        try:
            out.append(float(tok))
        except Exception:
            pass
    return out


def straglr_ref_len(v: cyvcf2.Variant) -> ty.Optional[float]:
    """
    Prefer RUL_REF (bp length in reference).
    Fallback: END-POS+1 if END exists.
    Last resort: len(REF).
    """
    rul_ref = v.INFO.get("RUL_REF")
    if rul_ref is not None:
        try:
            return float(rul_ref)
        except Exception:
            pass

    end = v.INFO.get("END")
    if end is not None:
        try:
            return float(int(end) - int(v.POS) + 1)
        except Exception:
            pass

    try:
        return float(len(v.REF))
    except Exception:
        return None


def straglr_allele_len_from_index(v: cyvcf2.Variant, idx: int) -> ty.Optional[float]:
    """
    idx=0 -> ref length (Lref)
    idx>0 -> ALT allele length from RB[idx-1]
    """
    if idx < 0:
        return None
    Lref = straglr_ref_len(v)
    if Lref is None:
        return None
    if idx == 0:
        return Lref

    rb_list = _parse_info_list(v.INFO.get("RB"))
    j = idx - 1
    if j < 0 or j >= len(rb_list):
        return None
    return float(rb_list[j])


def get_allele_lengths_straglr(v: cyvcf2.Variant) -> ty.Tuple[ty.Optional[float], ty.Optional[float], ty.Optional[float]]:
    """
    Returns (Lref, L1, L2) for the first sample based on GT and INFO.
    """
    Lref = straglr_ref_len(v)
    if Lref is None:
        return (None, None, None)

    try:
        g = v.genotypes[0]
        if len(g) < 2:
            return (Lref, None, None)
        a_idx, b_idx = int(g[0]), int(g[1])
    except Exception:
        return (Lref, None, None)

    L1 = straglr_allele_len_from_index(v, a_idx)
    L2 = straglr_allele_len_from_index(v, b_idx)
    return (Lref, L1, L2)


def length_features_from_lengths(Lref: float, L1: float, L2: float) -> ty.List[float]:
    """
    Feature vector in bp-space:
      [L1, L2, dL1, dL2] where dL = L - Lref
    """
    return [float(L1), float(L2), float(L1 - Lref), float(L2 - Lref)]


def _fmt_int_or_dot(x: ty.Optional[float]) -> str:
    if x is None:
        return "."
    try:
        return str(int(round(float(x))))
    except Exception:
        return "."


# -----------------------------
# Trio inheritance minimization (numeric lengths)
# -----------------------------
def parental_combos_len(
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
) -> ty.List[ty.Tuple[float, float, str]]:
    out = []
    for mi, di in product([0, 1], [0, 1]):
        out.append((mom[mi], dad[di], f"m{mi}d{di}"))
    return out


def _pick_tx_from_label_len(
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
    parent_ht: str,
) -> ty.Tuple[float, float]:
    # parent_ht like "m0d1"
    mi = int(parent_ht[1])
    di = int(parent_ht[3])
    return mom[mi], dad[di]


def min_inheritance_len_distance_numeric(
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
    kid: ty.Tuple[float, float],
    Lref: float,
    power: int = 1,
) -> ty.Tuple[float, str, ty.Tuple[float, float], str]:
    """
    Distance-only selection of best parental transmission and best kid ordering.

    Returns:
      (len_dist, parent_ht_label, kid_ht_pair, kid_ht_ordered)
    """
    kid_orders = [
        (kid[0], kid[1], "as_is"),
        (kid[1], kid[0], "flipped"),
    ]

    best = None  # (dist, lab, (k1,k2), order_tag)
    for k1, k2, order_tag in kid_orders:
        kid_vec = length_features_from_lengths(Lref, k1, k2)
        for m_a, d_a, lab in parental_combos_len(mom, dad):
            p_vec = length_features_from_lengths(Lref, m_a, d_a)
            dist = minkowski(p_vec, kid_vec, power=power)
            cand = (dist, lab, (k1, k2), order_tag)
            if best is None or cand[0] < best[0]:
                best = cand

    assert best is not None
    best_dist, best_lab, best_kid_pair, _ = best

    def best_dist_for(k1: float, k2: float) -> float:
        kid_vec = length_features_from_lengths(Lref, k1, k2)
        dmin = None
        for m_a, d_a, _lab in parental_combos_len(mom, dad):
            p_vec = length_features_from_lengths(Lref, m_a, d_a)
            d = minkowski(p_vec, kid_vec, power=power)
            dmin = d if dmin is None else min(dmin, d)
        return float(dmin) if dmin is not None else float("inf")

    d_as = best_dist_for(kid[0], kid[1])
    d_fl = best_dist_for(kid[1], kid[0])
    kid_ht_ordered = "T" if d_as != d_fl else "F"

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
      - merge_key="pos": (CHROM, POS)
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
    Straglr-compatible Mendelian check using INFO-derived allele lengths (RB/RUL_REF/END) + GT.
    Column NAMES are kept EXACTLY the same as your original output.
    For Straglr (no sequences), we WRITE BP LENGTHS into the existing *_seq columns:
      - mom_a1_seq/mom_a2_seq, dad_a1_seq/dad_a2_seq, kid_a1_seq/kid_a2_seq -> allele length (bp)
      - mom_tx_seq/dad_tx_seq -> transmitted allele length (bp)
    Everything else stays in the same columns.
    Consistency rule: len_dist == 0 -> consistent, else inconsistent.
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

    # KEEP COLUMN NAMES EXACTLY SAME
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

        # Extract lengths for parents
        mom_Lref, mom_L1, mom_L2 = get_allele_lengths_straglr(mom)
        dad_Lref, dad_L1, dad_L2 = get_allele_lengths_straglr(dad)

        if mom_Lref is None or dad_Lref is None:
            continue

        # Anchor reference length from mom (deterministic)
        Lref = float(mom_Lref)

        if mom_L1 is None or mom_L2 is None or dad_L1 is None or dad_L2 is None:
            continue

        mom_vec = length_features_from_lengths(Lref, float(mom_L1), float(mom_L2))
        dad_vec = length_features_from_lengths(Lref, float(dad_L1), float(dad_L2))

        mom_len = f"{int(mom_vec[0])},{int(mom_vec[1])}"
        mom_dlen = f"{int(mom_vec[2])},{int(mom_vec[3])}"
        dad_len = f"{int(dad_vec[0])},{int(dad_vec[1])}"
        dad_dlen = f"{int(dad_vec[2])},{int(dad_vec[3])}"

        # Put BP lengths into *_seq columns (keeps schema uniform)
        mom_a1_seq = _fmt_int_or_dot(mom_L1)
        mom_a2_seq = _fmt_int_or_dot(mom_L2)
        dad_a1_seq = _fmt_int_or_dot(dad_L1)
        dad_a2_seq = _fmt_int_or_dot(dad_L2)

        for i, kid in enumerate(kids):
            if kid is None:
                continue

            kid_Lref, kid_L1, kid_L2 = get_allele_lengths_straglr(kid)
            if kid_L1 is None or kid_L2 is None:
                continue

            kid_vec = length_features_from_lengths(Lref, float(kid_L1), float(kid_L2))
            kid_len = f"{int(kid_vec[0])},{int(kid_vec[1])}"
            kid_dlen = f"{int(kid_vec[2])},{int(kid_vec[3])}"

            len_dist, parent_ht, kid_ht_pair, kid_ht_ordered = min_inheritance_len_distance_numeric(
                (float(mom_L1), float(mom_L2)),
                (float(dad_L1), float(dad_L2)),
                (float(kid_L1), float(kid_L2)),
                Lref=Lref,
                power=power,
            )

            mom_tx_len, dad_tx_len = _pick_tx_from_label_len(
                (float(mom_L1), float(mom_L2)),
                (float(dad_L1), float(dad_L2)),
                parent_ht,
            )

            p_vec = length_features_from_lengths(Lref, mom_tx_len, dad_tx_len)
            k_vec = length_features_from_lengths(Lref, kid_ht_pair[0], kid_ht_pair[1])

            absdiff = [abs(p_vec[j] - k_vec[j]) for j in range(4)]
            dist_sum = minkowski_sum_only(p_vec, k_vec, power=power)

            mendelian = "T" if len_dist == 0 else "F"
            total_scored += 1
            if mendelian == "T":
                total_consistent += 1
            len_dists_all.append(float(len_dist))

            # kid_ht stays as "a,b" (bp lengths)
            kid_ht = f"{int(kid_ht_pair[0])},{int(kid_ht_pair[1])}"

            # kid_ID from VCF ID column
            kid_ID = kid.ID if kid.ID is not None and kid.ID != "" else "."

            # Put BP lengths into kid *_seq columns and transmitted columns
            kid_a1_seq = _fmt_int_or_dot(kid_L1)
            kid_a2_seq = _fmt_int_or_dot(kid_L2)
            mom_tx_seq = _fmt_int_or_dot(mom_tx_len)
            dad_tx_seq = _fmt_int_or_dot(dad_tx_len)

            line = (
                f"{chrom}\t{pos}\t{kid_ids[i]}\t{kid_ID}\t"
                f"{fmt_gt(mom)}\t{fmt_gt(dad)}\t{fmt_gt(kid)}\t"
                f"{mom.REF}\t{_alts_str(mom)}\t"
                f"{dad.REF}\t{_alts_str(dad)}\t"
                f"{kid.REF}\t{_alts_str(kid)}\t"
                f"{mom_a1_seq}\t{mom_a2_seq}\t"
                f"{dad_a1_seq}\t{dad_a2_seq}\t"
                f"{kid_a1_seq}\t{kid_a2_seq}\t"
                f"{int(Lref)}\t"
                f"{mom_len}\t{mom_dlen}\t"
                f"{dad_len}\t{dad_dlen}\t"
                f"{kid_len}\t{kid_dlen}\t"
                f"{parent_ht}\t{mom_tx_seq}\t{dad_tx_seq}\t"
                f"{kid_ht}\t{kid_ht_ordered}\t"
                f"{int(absdiff[0])}\t{int(absdiff[1])}\t{int(absdiff[2])}\t{int(absdiff[3])}\t"
                f"{dist_sum}\t{len_dist}\t{mendelian}"
            )
            print(line, file=fh)

    fh.close()

    pct = (total_consistent / total_scored) * 100.0 if total_scored > 0 else 0.0
    print(
        f"[summary] scored_loci={total_scored} consistent={total_consistent} pct={pct:.2f}%",
        file=sys.stderr,
    )

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
        description=(
            "Straglr-compatible Mendelian consistency using INFO-derived allele lengths (RB/RUL_REF/END) + GT. "
            "Keeps ORIGINAL column names; writes BP lengths into *_seq and *_tx_seq value fields."
        )
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
        help="How to merge loci across VCFs: 'pos' (CHROM+POS, default) or 'id' (ID only).",
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
