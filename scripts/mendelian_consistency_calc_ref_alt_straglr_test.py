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
    if merge_key == "id":
        return v.ID if v.ID is not None and v.ID != "" else "."
    return (chrom_key(v.CHROM), int(v.POS))


# -----------------------------
# Distance functions
# -----------------------------
def minkowski(a: ty.List[float], b: ty.List[float], power: int = 1) -> float:
    return (sum(abs(a[i] - b[i]) ** power for i in range(len(a))) ** (1.0 / power))


def minkowski_sum_only(a: ty.List[float], b: ty.List[float], power: int = 1) -> float:
    return sum(abs(a[i] - b[i]) ** power for i in range(len(a)))


# -----------------------------
# GT formatting
# -----------------------------
def fmt_gt(v: cyvcf2.Variant) -> str:
    """
    Return GT string. Handles both diploid (0/1) and haploid (0, 1) Straglr GTs.
    """
    try:
        g = v.genotypes[0]
        if len(g) < 1:
            return "."
        a = g[0]
        # Haploid: cyvcf2 fills second allele as -1 for haploid calls
        b = g[1] if len(g) >= 2 else -1
        if a is None or a < 0:
            return "."
        if b is None or b < 0:
            return f"{a}"
        phased = bool(g[2]) if len(g) >= 3 else False
        sep = "|" if phased else "/"
        return f"{a}{sep}{b}"
    except Exception:
        return "."


def _alts_str(v: cyvcf2.Variant) -> str:
    alts = list(v.ALT) if v.ALT is not None else []
    return ",".join(alts) if alts else "."


def _fmt_int_or_dot(x: ty.Optional[float]) -> str:
    if x is None:
        return "."
    try:
        return str(int(round(float(x))))
    except Exception:
        return "."


# -----------------------------
# FIX 1: Straglr length extraction
# -----------------------------
def _parse_info_list(val) -> ty.List[float]:
    """
    Parse cyvcf2 INFO values into list[float].
    Handles None, int/float, string "a,b,c", tuple/list.
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
    # Priority 1: END - POS + 1 (correct locus span)
    end = v.INFO.get("END")
    if end is not None:
        try:
            return float(int(end) - int(v.POS) + 1)
        except Exception:
            pass

    # Priority 2: len(REF) sequence
    try:
        return float(len(v.REF))
    except Exception:
        return None


def straglr_allele_len_from_index(
    v: cyvcf2.Variant,
    idx: int,
    Lref: float,
) -> ty.Optional[float]:
    """
    Resolve a single allele length from its GT index.

    Straglr encodes allele lengths as follows:
      - REF allele (idx=0): length = END - POS + 1  (the reference locus span = Lref)
      - ALT allele (idx>0): length = RB[idx-1]      (bp count from INFO/RB field)
      - Missing   (idx<0):  -> None                 (haploid/missing handled upstream)

    RB is a comma-separated list of ALT allele lengths in bp, one per ALT:
      e.g.  GT=0/1  RB=273        -> ALT allele = 273 bp,  REF allele = Lref
            GT=1/2  RB=310,684    -> ALT1 = 310 bp, ALT2 = 684 bp
            GT=0/0  (no RB)       -> both alleles = Lref
            GT=1    (haploid alt) -> RB[0] duplicated (handled in get_allele_lengths_straglr)

    Per-GT mapping summary:
      GT=0/0  -> (Lref,    Lref   )   both from END-POS+1
      GT=0/1  -> (Lref,    RB[0]  )   ref from END-POS+1, alt from RB
      GT=1/1  -> (RB[0],   RB[0]  )   both from RB[0]
      GT=1/2  -> (RB[0],   RB[1]  )   each alt from its own RB entry
      GT=0    -> (Lref,    Lref   )   haploid ref -> treated as 0/0
      GT=1    -> (RB[0],   RB[0]  )   haploid alt -> treated as 1/1
    """
    if idx < 0:
        return None
    if idx == 0:
        return Lref
    rb_list = _parse_info_list(v.INFO.get("RB"))
    j = idx - 1
    if j < 0 or j >= len(rb_list):
        return None
    return float(rb_list[j])


def get_allele_lengths_straglr(
    v: cyvcf2.Variant,
    anchor_Lref: ty.Optional[float] = None,
) -> ty.Tuple[ty.Optional[float], ty.Optional[float], ty.Optional[float]]:
    """
    Returns (Lref, L1, L2) for the first sample using GT + INFO fields.

    Length sources per GT class
    ---------------------------
    GT=0/0  no RB field expected
            L1 = END-POS+1  (Lref)
            L2 = END-POS+1  (Lref)
            Example: END=648889, POS=648806 -> Lref=84, L1=84, L2=84

    GT=0/1  RB has 1 value (the single ALT allele)
            L1 = END-POS+1  (Lref, REF allele)
            L2 = RB[0]      (ALT allele bp length)
            Example: END=10800, POS=10627, RB=273 -> Lref=174, L1=174, L2=273

    GT=1/1  RB has 1 value (both alleles same ALT)
            L1 = RB[0]
            L2 = RB[0]

    GT=1/2  RB has 2 values (one per ALT allele)
            L1 = RB[0]
            L2 = RB[1]
            Example: END=181463, POS=181146, RB=310,684 -> L1=310, L2=684

    GT=0    (haploid ref — Straglr-specific, no slash)
            Treated as 0/0: L1=Lref, L2=Lref

    GT=1    (haploid alt — Straglr-specific, no slash)
            Treated as 1/1: L1=RB[0], L2=RB[0]

    anchor_Lref: if provided (from mom), overrides this variant's own Lref
                 so all trio members use a consistent reference span.
    """
    Lref = straglr_ref_len(v)
    if Lref is None:
        return (None, None, None)

    # Use anchor Lref from mom if provided (keeps trio comparison consistent)
    effective_Lref = anchor_Lref if anchor_Lref is not None else Lref

    try:
        g = v.genotypes[0]
        if len(g) < 1:
            return (effective_Lref, None, None)

        a_idx = int(g[0])
        # FIX 2: detect haploid (second allele = -1 or missing)
        b_raw = g[1] if len(g) >= 2 else -1
        b_idx = int(b_raw) if (b_raw is not None and int(b_raw) >= 0) else None

        if b_idx is None:
            # Haploid call: duplicate the single allele -> treat as homozygous
            L = straglr_allele_len_from_index(v, a_idx, effective_Lref)
            return (effective_Lref, L, L)

        L1 = straglr_allele_len_from_index(v, a_idx, effective_Lref)
        L2 = straglr_allele_len_from_index(v, b_idx, effective_Lref)
        return (effective_Lref, L1, L2)

    except Exception:
        return (effective_Lref, None, None)


# -----------------------------
# FIX 3: Simplified feature vector [L1, L2] only
# -----------------------------
def length_features(L1: float, L2: float) -> ty.List[float]:
    """
    FIX 3: Use [L1, L2] only — same as fixed general script.

    OLD: [L1, L2, dL1, dL2] was redundant (dL = L - Lref, constant offset)
         and caused float noise in dL to push exact matches to non-zero distance.
    FIXED: [L1, L2] is the minimal sufficient feature vector for Mendelian
           consistency — the only question is whether kid lengths match a
           parental combination.
    """
    return [float(L1), float(L2)]


# -----------------------------
# Trio inheritance minimization
# -----------------------------
def parental_combos(
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
) -> ty.List[ty.Tuple[float, float, str]]:
    out = []
    for mi, di in product([0, 1], [0, 1]):
        out.append((mom[mi], dad[di], f"m{mi}d{di}"))
    return out


def _pick_tx_from_label(
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
    parent_ht: str,
) -> ty.Tuple[float, float]:
    mi = int(parent_ht[1])
    di = int(parent_ht[3])
    return mom[mi], dad[di]


def best_dist_for_ordering(
    k1: float,
    k2: float,
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
    power: int,
) -> float:
    """Minimum distance over all parental combos for a given kid ordering."""
    kid_vec = length_features(k1, k2)
    dmin = None
    for m_a, d_a, _lab in parental_combos(mom, dad):
        p_vec = length_features(m_a, d_a)
        d = minkowski(p_vec, kid_vec, power=power)
        dmin = d if dmin is None else min(dmin, d)
    return float(dmin) if dmin is not None else float("inf")


def min_inheritance_distance(
    mom: ty.Tuple[float, float],
    dad: ty.Tuple[float, float],
    kid: ty.Tuple[float, float],
    power: int = 1,
) -> ty.Tuple[float, str, ty.Tuple[float, float], str]:
    """
    Find the minimum-distance parental transmission and best kid allele ordering.

    Returns:
      (best_dist, parent_ht_label, kid_ht_pair, kid_ht_ordered)
    """
    kid_orders = [
        (kid[0], kid[1], "as_is"),
        (kid[1], kid[0], "flipped"),
    ]

    best = None
    for k1, k2, order_tag in kid_orders:
        kid_vec = length_features(k1, k2)
        for m_a, d_a, lab in parental_combos(mom, dad):
            p_vec = length_features(m_a, d_a)
            dist = minkowski(p_vec, kid_vec, power=power)
            cand = (dist, lab, (k1, k2), order_tag)
            if best is None or cand[0] < best[0]:
                best = cand

    assert best is not None
    best_dist, best_lab, best_kid_pair, _ = best

    # Determine if kid allele phase is resolvable
    d_as = best_dist_for_ordering(kid[0], kid[1], mom, dad, power)
    d_fl = best_dist_for_ordering(kid[1], kid[0], mom, dad, power)
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
    output_prefix: str = "straglr-mc-",
    exclude_chroms: ty.List[str] = None,
    merge_key: str = "pos",
    consistency_threshold: float = 0.5,
):
    """
    Straglr-specific Mendelian consistency script.

    Key fixes vs previous version:
      FIX 1: RUL_REF is motif size NOT locus length. Lref now = END-POS+1.
      FIX 2: Haploid GTs (0, 1) are now handled — treated as homozygous.
      FIX 3: Feature vector simplified to [L1, L2] — removes float noise from dL.
      FIX 4: Float-safe consistency threshold (default 0.5 bp) instead of == 0.
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
        "Lref\t"
        "mom_L1\tmom_L2\t"
        "dad_L1\tdad_L2\t"
        "kid_L1\tkid_L2\t"
        "parent_ht\tmom_tx_len\tdad_tx_len\t"
        "kid_ht\tkid_ht_ordered\t"
        "absdiff_L1\tabsdiff_L2\t"
        "distance_power_sum\tlen_dist\tmendelian\t"
        "mom_haploid\tdad_haploid\tkid_haploid"
    )
    print(header, file=fh)

    n_kids = len(kid_vcfs)
    total_scored = [0] * n_kids
    total_consistent = [0] * n_kids
    total_skipped_no_len = [0] * n_kids
    total_haploid_kid = [0] * n_kids

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

        # FIX 1: Lref from END-POS+1, anchored to mom
        mom_Lref, mom_L1, mom_L2 = get_allele_lengths_straglr(mom)
        if mom_Lref is None or mom_L1 is None or mom_L2 is None:
            continue

        Lref = float(mom_Lref)

        # Pass anchor Lref to dad and kid for consistent dL computation
        dad_Lref, dad_L1, dad_L2 = get_allele_lengths_straglr(dad, anchor_Lref=Lref)
        if dad_L1 is None or dad_L2 is None:
            continue

        # Detect haploid calls for reporting
        def is_haploid(v: cyvcf2.Variant) -> bool:
            try:
                g = v.genotypes[0]
                b = g[1] if len(g) >= 2 else -1
                return b is None or int(b) < 0
            except Exception:
                return False

        mom_haploid = is_haploid(mom)
        dad_haploid = is_haploid(dad)

        for i, kid in enumerate(kids):
            if kid is None:
                continue

            kid_haploid = is_haploid(kid)

            # FIX 2: pass anchor Lref so kid uses same reference span as mom
            _, kid_L1, kid_L2 = get_allele_lengths_straglr(kid, anchor_Lref=Lref)

            if kid_L1 is None or kid_L2 is None:
                total_skipped_no_len[i] += 1
                continue

            if kid_haploid:
                total_haploid_kid[i] += 1

            len_dist, parent_ht, kid_ht_pair, kid_ht_ordered = min_inheritance_distance(
                (float(mom_L1), float(mom_L2)),
                (float(dad_L1), float(dad_L2)),
                (float(kid_L1), float(kid_L2)),
                power=power,
            )

            mom_tx_len, dad_tx_len = _pick_tx_from_label(
                (float(mom_L1), float(mom_L2)),
                (float(dad_L1), float(dad_L2)),
                parent_ht,
            )

            p_vec = length_features(mom_tx_len, dad_tx_len)
            k_vec = length_features(kid_ht_pair[0], kid_ht_pair[1])

            absdiff = [abs(p_vec[j] - k_vec[j]) for j in range(2)]
            dist_sum = minkowski_sum_only(p_vec, k_vec, power=power)

            # FIX 4: float-safe threshold (default 0.5 bp)
            # This prevents floating point noise from wrongly calling inconsistency
            # and also correctly populates the ONE-OFF bucket downstream
            mendelian = "T" if len_dist <= consistency_threshold else "F"

            total_scored[i] += 1
            if mendelian == "T":
                total_consistent[i] += 1
            len_dists_all.append(float(len_dist))

            kid_ht = f"{int(round(kid_ht_pair[0]))},{int(round(kid_ht_pair[1]))}"
            kid_ID = kid.ID if kid.ID is not None and kid.ID != "" else "."

            line = (
                f"{chrom}\t{pos}\t{kid_ids[i]}\t{kid_ID}\t"
                f"{fmt_gt(mom)}\t{fmt_gt(dad)}\t{fmt_gt(kid)}\t"
                f"{mom.REF}\t{_alts_str(mom)}\t"
                f"{dad.REF}\t{_alts_str(dad)}\t"
                f"{kid.REF}\t{_alts_str(kid)}\t"
                f"{int(Lref)}\t"
                f"{_fmt_int_or_dot(mom_L1)}\t{_fmt_int_or_dot(mom_L2)}\t"
                f"{_fmt_int_or_dot(dad_L1)}\t{_fmt_int_or_dot(dad_L2)}\t"
                f"{_fmt_int_or_dot(kid_L1)}\t{_fmt_int_or_dot(kid_L2)}\t"
                f"{parent_ht}\t{_fmt_int_or_dot(mom_tx_len)}\t{_fmt_int_or_dot(dad_tx_len)}\t"
                f"{kid_ht}\t{kid_ht_ordered}\t"
                f"{int(round(absdiff[0]))}\t{int(round(absdiff[1]))}\t"
                f"{dist_sum:.4f}\t{len_dist:.4f}\t{mendelian}\t"
                f"{int(mom_haploid)}\t{int(dad_haploid)}\t{int(kid_haploid)}"
            )
            print(line, file=fh)

    fh.close()

    # Per-kid summary
    print("[summary]", file=sys.stderr)
    for i, kid_id in enumerate(kid_ids):
        scored = total_scored[i]
        consistent = total_consistent[i]
        skipped = total_skipped_no_len[i]
        haploid = total_haploid_kid[i]
        pct = (consistent / scored * 100.0) if scored > 0 else 0.0
        print(
            f"  kid={kid_id} scored={scored} consistent={consistent} "
            f"pct={pct:.2f}% skipped_no_len={skipped} haploid_treated_as_hom={haploid}",
            file=sys.stderr,
        )

    if len(len_dists_all) > 0:
        mname = "manhattan" if power == 1 else ("euclidean" if power == 2 else f"minkowski(p={power})")
        fig = px.histogram(x=len_dists_all)
        fig.update_layout(
            xaxis_title=f"{mname} distance [L1,L2] Straglr (min transmission)",
            yaxis_title="count",
        )
        fig.write_html(out_html)
        print(f"wrote {out_html}", file=sys.stderr)

    print(f"wrote {out_txt}", file=sys.stderr)


def parse_args(argv=None):
    import argparse

    p = argparse.ArgumentParser(
        description="Straglr Mendelian consistency using INFO-derived lengths (END/RB) + GT."
    )
    p.add_argument("--mom", required=True, type=pathlib.Path)
    p.add_argument("--dad", required=True, type=pathlib.Path)
    p.add_argument(
        "--kids",
        required=True,
        action="append",
        type=pathlib.Path,
        help="Repeat --kids for multiple kids.",
    )
    p.add_argument("--power", type=int, default=1, help="1=Manhattan, 2=Euclidean (default: 1).")
    p.add_argument("--out-prefix", default="straglr-mc-", help="Output prefix.")
    p.add_argument(
        "--exclude-chroms",
        nargs="*",
        default=["chrX", "chrY"],
    )
    p.add_argument(
        "--merge-key",
        choices=["pos", "id"],
        default="pos",
    )
    p.add_argument(
        "--consistency-threshold",
        type=float,
        default=0.5,
        help=(
            "Max distance to call consistent (default: 0.5 bp). "
            "Prevents float precision errors from calling exact matches inconsistent. "
            "Use 1.0 to also count 1bp-off as consistent."
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