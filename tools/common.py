# -*- coding: utf-8 -*-
"""Shared oligo-splitting logic used by Silent Library Designer and Oligo Assembly Designer."""

COMPLEMENTARY_TABLE = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "D": "H",
    "H": "D", "V": "B", "N": "N",
}

MIN_OLIGO_LEN = 50
MAX_OLIGO_LEN = 90
NOMINAL_OVERLAP_FOR_PLANNING = 20  # auto mode no longer takes overlap from the user - this is
# only a seed for the initial oligo-count/spacing plan; the actual overlap at each junction is
# freely re-optimized for Tm afterward and reported back per candidate


def reverse_complement(sequence):
    return "".join(COMPLEMENTARY_TABLE.get(base, base) for base in sequence)[::-1]


def create_oligos(sequence, oligo_length, overlap):
    """Tile `sequence` into overlapping oligos of `oligo_length` with `overlap` shared bases.
    The final oligo may be shorter than `oligo_length`."""
    step = oligo_length - overlap
    oligos = []
    start = 0
    n = len(sequence)
    while start < n:
        end = min(start + oligo_length, n)
        oligos.append((start + 1, end, sequence[start:end]))
        if end == n:
            break
        start += step
    return oligos


def plan_oligo_split(length, overlap):
    """Find the smallest even oligo count whose length fits MIN_OLIGO_LEN..MAX_OLIGO_LEN
    while keeping ~`overlap` shared bases between consecutive oligos.

    Returns (num_oligos, oligo_length) or None if no such split exists.
    """
    if length <= 0 or overlap < 0 or overlap >= MAX_OLIGO_LEN or length <= overlap:
        return None

    num_oligos = 2
    while True:
        avg_length = overlap + (length - overlap) / num_oligos
        if avg_length <= MAX_OLIGO_LEN:
            break
        num_oligos += 2

    if avg_length < MIN_OLIGO_LEN:
        return None

    oligo_length = max(MIN_OLIGO_LEN, min(MAX_OLIGO_LEN, round(avg_length)))
    return num_oligos, oligo_length


def create_balanced_oligos(sequence, num_oligos, oligo_length):
    """Tile `sequence` into `num_oligos` oligos of `oligo_length`, spaced evenly
    so the (slightly varying) overlap is spread out instead of dumped on the last piece."""
    n = len(sequence)
    if num_oligos <= 1 or n <= oligo_length:
        return [(1, n, sequence)]

    oligos = []
    for i in range(num_oligos):
        start = round(i * (n - oligo_length) / (num_oligos - 1))
        end = start + oligo_length
        oligos.append((start + 1, end, sequence[start:end]))
    return oligos




# --- PCA-style (annealing) oligo assembly: alternating F/R strand, Tm-driven overlaps ---
# Used by Oligo Assembly Designer only (Silent Library's IUPAC-degenerate sequence has no
# single well-defined Tm, so it keeps the fixed-length tiling above).
#
# Overlap design parameters follow standard PCA (polymerase cycling assembly) practice:
# 18-25bp overlaps, 40-60% GC, Tm 55-57C (52C hard floor), a single G/C 3' clamp preferred
# (not a run of them), high-T oligos flagged (T-richness promotes mispairing).

MIN_OVERLAP = 18
MAX_OVERLAP = 25
EXTENDED_MAX_OVERLAP = 90  # fallback ceiling when 18-25bp can't reach the Tm floor (e.g. AT-rich regions) - capped at the max synthesizable oligo length (90mer)
JUNCTION_SHIFT_TOLERANCE = 10  # how far a junction may slide off its nominal midpoint
TARGET_TM = 56.0  # user-adjustable "roughly" via the Oligo Assembly Designer form; default reproduces the 55-57C sweet spot / 52C floor protocol spec
TM_SWEET_SPREAD = 1.0  # sweet spot = target +/- this
TM_FLOOR_MARGIN = 4.0  # hard floor = target - this
MAX_TM_DEVIATION = 4.0  # a design is unusable if any junction's Tm strays more than this from
# the design's own average junction Tm - PCA runs at one annealing temperature for every
# junction at once, so a single outlier can fail to anneal (or over-anneal) relative to the rest
TM_SEARCH_OFFSETS = [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5]  # instead of forcing the user's exact
# target Tm, search nearby values too and surface whichever actually satisfies the most rules
TM_SCAN_MIN = 45.0  # no Tm input from the user - scan this whole practical PCR/PCA range
TM_SCAN_MAX = 65.0
TM_SCAN_STEP = 2.0
MAX_CANDIDATES_SHOWN = 5  # cap on how many rows the candidate-comparison table displays
PRICE_PER_NT_KRW = 250  # rough oligo synthesis cost estimate, KRW per nt
TM_COMFORT_MIN = 50.0  # below this average junction Tm, overlap annealing gets progressively
# less robust in practice (low fragment concentration, low-GC overlaps, secondary structure,
# poor transformation efficiency all compound) - penalize it even though such a design may
# still be "usable" by the strict consistency check
TM_LOW_PENALTY_WEIGHT = 3.0  # score points per degree C below TM_COMFORT_MIN
GC_MIN = 0.40
GC_MAX = 0.60
HIGH_T_FRACTION = 0.40
HOMOPOLYMER_RUN = 5  # >= this many identical bases in a row = slippage/mispriming risk
GC_BASES = "GC"

# 3' terminus is the single most important region for PCR/PCA success and specificity:
# check its GC content, and scan for 3'-anchored dimers (cross-dimer with other oligos,
# self-dimer/hairpin, or mispriming elsewhere in the template) using real nearest-neighbor
# thermodynamics rather than a simple GC% proxy.
TERMINAL_WINDOW = 6  # the critical last N bases of the 3' end (5-8nt rule of thumb)
TERMINAL_GC_MIN = 2  # want at least this many G/C in the terminal window
DIMER_MIN_RUN = 4  # shorter runs are not thermodynamically meaningful to flag
DIMER_MAX_RUN = 12  # only the last ~12nt actually matters for a 3'-anchored dimer
DIMER_DG_CAUTION = -6.0  # kcal/mol: more stable (more negative) than this -> caution
DIMER_DG_REDESIGN = -9.0  # kcal/mol: more stable than this -> redesign this primer

# SantaLucia (1998) unified nearest-neighbor free energies at 37C, kcal/mol (Watson-Crick).
NN_DG37 = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}
DUPLEX_INIT_DG37 = 1.96  # approximate duplex-initiation penalty (terminal-independent)


def _has_homopolymer_run(seq, run_len=HOMOPOLYMER_RUN):
    run_base, run_count = None, 0
    for base in seq:
        if base == run_base:
            run_count += 1
        else:
            run_base, run_count = base, 1
        if run_count >= run_len:
            return True
    return False


def _occurs_elsewhere(sequence, seg, at_pos):
    """True if `seg` (or its reverse complement) appears anywhere in `sequence` other than
    at `at_pos` - i.e. this exact junction sequence is ambiguous and could let the oligo
    mis-anneal to the wrong partner instead of its intended neighbor."""
    if len(seg) < 6:
        return False
    rc = reverse_complement(seg)
    for needle in {seg, rc}:
        idx = sequence.find(needle)
        while idx != -1:
            if idx != at_pos:
                return True
            idx = sequence.find(needle, idx + 1)
    return False


def calc_duplex_dg(seq):
    """Approximate Gibbs free energy (kcal/mol, 37C) of the fully-matched duplex `seq` would
    form with its own perfect complement, via nearest-neighbor stacking (SantaLucia 1998
    unified parameters) plus a duplex-initiation penalty. More negative = more stable."""
    if len(seq) < 2:
        return 0.0
    dg = DUPLEX_INIT_DG37
    for i in range(len(seq) - 1):
        dg += NN_DG37.get(seq[i:i + 2], -1.3)
    return round(dg, 2)


def find_3prime_dimer(order_sequence, pool, min_run=DIMER_MIN_RUN, max_run=DIMER_MAX_RUN):
    """Scan for the longest run, anchored at the 3' terminus of `order_sequence`, whose reverse
    complement appears in any sequence of `pool` - i.e. a stretch capable of annealing (and
    being extended by polymerase, since it includes the 3'-most base) to something it
    shouldn't: another oligo (primer-dimer), itself (hairpin/self-dimer), or elsewhere in the
    template (mispriming). Returns (run_length, dg, source_index) for the longest such run, or
    None if nothing at or above min_run was found. `pool` is a list of (label, sequence)."""
    three_prime = order_sequence[-max_run:]
    for k in range(min(max_run, len(three_prime)), min_run - 1, -1):
        suffix = three_prime[-k:]
        rc = reverse_complement(suffix)
        for label, pool_seq in pool:
            if rc in pool_seq:
                return k, calc_duplex_dg(suffix), label
    return None


def dimer_flag(dg):
    if dg is None:
        return "ok"
    if dg <= DIMER_DG_REDESIGN:
        return "redesign"
    if dg <= DIMER_DG_CAUTION:
        return "caution"
    return "ok"


def _find_3prime_template_dimer(order_sequence, sequence, own_start0, own_end0,
                                 min_run=DIMER_MIN_RUN, max_run=DIMER_MAX_RUN):
    """Like find_3prime_dimer, but scanned against the full template while ignoring the
    oligo's own intended binding span (so its designed match doesn't count as 'elsewhere')."""
    three_prime = order_sequence[-max_run:]
    for k in range(min(max_run, len(three_prime)), min_run - 1, -1):
        suffix = three_prime[-k:]
        rc = reverse_complement(suffix)
        idx = sequence.find(rc)
        while idx != -1:
            if idx + k <= own_start0 or idx >= own_end0:
                return k, calc_duplex_dg(suffix), "template 다른 위치"
            idx = sequence.find(rc, idx + 1)
    return None


def _dimer_check_for_row(i, rows, sequence):
    """Worst (most stable / most negative dG) 3'-anchored dimer found for row i, checked
    against every OTHER, NON-ADJACENT oligo (primer-dimer), itself (hairpin/self-dimer), and
    the template elsewhere (mispriming) - or None if nothing at or above DIMER_MIN_RUN was
    found. Immediate neighbors (i-1, i+1) are excluded from the cross-dimer pool: they are
    *designed* to be complementary at their shared junction, so a match there is PCA working
    correctly, not a primer-dimer defect."""
    row = rows[i]
    order_seq = row["order_sequence"]
    other_pool = [
        (f"올리고 #{r['index']}", r["order_sequence"])
        for j, r in enumerate(rows)
        if j != i and abs(j - i) > 1
    ]
    self_pool = [("자기 자신(헤어핀/자가이량체)", order_seq)]

    findings = []
    found = find_3prime_dimer(order_seq, other_pool)
    if found:
        findings.append(found)
    found = find_3prime_dimer(order_seq, self_pool)
    if found:
        findings.append(found)
    found = _find_3prime_template_dimer(order_seq, sequence, row["start"] - 1, row["end"])
    if found:
        findings.append(found)

    return min(findings, key=lambda f: f[1]) if findings else None


def add_primer_qc(rows, sequence):
    """Attach 3'-terminus quality fields to each row in place: dimer_run/dimer_dg/
    dimer_source/dimer_flag (cross-dimer, self-dimer/hairpin, or template mispriming - see
    find_3prime_dimer) and terminal_gc_count/terminal_gc_ok (GC content of the last
    TERMINAL_WINDOW bases, the single most important stretch for PCR specificity)."""
    for i, row in enumerate(rows):
        found = _dimer_check_for_row(i, rows, sequence)
        if found:
            run_len, dg, source = found
            row["dimer_run"], row["dimer_dg"], row["dimer_source"] = run_len, dg, source
        else:
            row["dimer_run"] = row["dimer_dg"] = row["dimer_source"] = None
        row["dimer_flag"] = dimer_flag(row["dimer_dg"])

        term_window = row["order_sequence"][-TERMINAL_WINDOW:]
        row["terminal_gc_count"] = sum(1 for b in term_window if b in GC_BASES)
        row["terminal_gc_ok"] = row["terminal_gc_count"] >= TERMINAL_GC_MIN
    return rows


def calc_tm(seq):
    """Rough oligo melting-temperature estimate, no salt/Mg correction:
    Wallace rule below 14nt, the %GC method (IDT/Promega convention) at or above it."""
    n = len(seq)
    if n == 0:
        return 0.0
    gc = sum(1 for base in seq if base in "GC")
    at = n - gc
    if n < 14:
        return float(2 * at + 4 * gc)
    return 64.9 + 41 * (gc - 16.4) / n


def calc_gc_fraction(seq):
    if not seq:
        return 0.0
    return sum(1 for base in seq if base in "GC") / len(seq)


def tm_bounds(target_tm):
    """(sweet_min, sweet_max, floor) around a user-chosen target Tm."""
    return target_tm - TM_SWEET_SPREAD, target_tm + TM_SWEET_SPREAD, target_tm - TM_FLOOR_MARGIN


def _junction_score(sequence, a, b, offset=0, target_tm=TARGET_TM, controls_termini=True):
    """Lower is better. Primary: land the Tm in the target +/-1C sweet spot (hard-floored at
    target-4C). Secondary: keep overlap GC% within 40-60%. Tie-break: since oligos alternate
    strand, this junction sets a real 3' end on BOTH flanking oligos only when the left one is
    forward (`controls_termini`) - the forward oligo's 3' end is sequence[b-1], and the reverse
    oligo's 3' end is complement(sequence[a]) which is G/C exactly when sequence[a] is. Prefer
    a single G/C there (a "GC clamp") but not two in a row (over-stabilizes). On the other
    alternation (left=reverse/right=forward) this junction only sets 5' ends, so the clamp
    check is skipped. Also steers the junction away from homopolymer runs / sequence repeated
    elsewhere in the gene, so any such risky stretch lands in the middle of an oligo instead of
    an annealing boundary."""
    sweet_min, sweet_max, floor = tm_bounds(target_tm)
    seg = sequence[a:b]
    tm = calc_tm(seg)
    gc_frac = calc_gc_fraction(seg)

    if tm < sweet_min:
        tm_penalty = sweet_min - tm
    elif tm > sweet_max:
        tm_penalty = tm - sweet_max
    else:
        tm_penalty = 0.0
    floor_penalty = 0.0 if tm >= floor else (floor - tm) * 50

    if gc_frac < GC_MIN:
        gc_penalty = (GC_MIN - gc_frac) * 10
    elif gc_frac > GC_MAX:
        gc_penalty = (gc_frac - GC_MAX) * 10
    else:
        gc_penalty = 0.0

    terminus_penalty = 0.0
    if controls_termini:
        if sequence[b - 1] not in GC_BASES:
            terminus_penalty += 0.1
        elif b - 2 >= 0 and sequence[b - 2] in GC_BASES:
            terminus_penalty += 0.05
        if sequence[a] not in GC_BASES:
            terminus_penalty += 0.1
        elif a + 1 < len(sequence) and sequence[a + 1] in GC_BASES:
            terminus_penalty += 0.05

    risk_penalty = 0.0
    if _has_homopolymer_run(seg):
        risk_penalty += 2.0
    if _occurs_elsewhere(sequence, seg, a):
        risk_penalty += 2.0

    # Oligo/overlap length is intentionally NOT penalized here - it's free to grow as far as
    # MIN_OLIGO_LEN..MAX_OLIGO_LEN allows if that's what it takes to hit the Tm target.
    drift_penalty = abs(offset) * 0.15  # pull back toward the balanced nominal position, so
    # drift only happens when it meaningfully helps Tm/GC/risk - keeps later junctions from
    # being squeezed into infeasibility by an earlier junction drifting too far

    # Tm match is the top priority (explicit user preference), weighted well above the
    # GC%/GC-clamp/repeat-risk tie-breakers.
    return (
        tm_penalty * 8 + floor_penalty
        + gc_penalty + terminus_penalty + risk_penalty + drift_penalty
    ), tm


def _find_junction_cut(sequence, nominal_start, nominal_end,
                        min_overlap=MIN_OVERLAP, max_overlap=EXTENDED_MAX_OVERLAP, target_tm=TARGET_TM,
                        left_start=None, right_end_fixed=None, controls_termini=True):
    """Search overlap lengths min_overlap..max_overlap (by default the full MIN_OVERLAP..
    EXTENDED_MAX_OVERLAP range - overlap/oligo length is unconstrained, so there's no reason to
    prefer a shorter one), and positions within +/-JUNCTION_SHIFT_TOLERANCE of the nominal
    junction midpoint, for the (start, end, length, tm) that best satisfies the PCA design
    targets around `target_tm`. A homopolymer run or a sequence repeated elsewhere in the gene
    is penalized so the search prefers to slide the junction into a cleaner spot rather than
    sit on top of it.

    `left_start` (the already-fixed start of the oligo ending at this junction) and
    `right_end_fixed` (set only for the last junction, whose right oligo end is the sequence
    end) keep both flanking oligos within MIN_OLIGO_LEN..MAX_OLIGO_LEN - the only length
    constraint that still applies is the 90nt physical synthesis cap."""
    mid = round((nominal_start + nominal_end) / 2)

    def search(lo, hi):
        best = None
        for offset in range(-JUNCTION_SHIFT_TOLERANCE, JUNCTION_SHIFT_TOLERANCE + 1):
            center = mid + offset
            for k in range(lo, hi + 1):
                a = center - k // 2
                b = a + k
                if a < 0 or b > len(sequence):
                    continue
                if left_start is not None and not (MIN_OLIGO_LEN <= b - left_start <= MAX_OLIGO_LEN):
                    continue
                if right_end_fixed is not None and not (MIN_OLIGO_LEN <= right_end_fixed - a <= MAX_OLIGO_LEN):
                    continue
                score, tm = _junction_score(sequence, a, b, offset, target_tm, controls_termini)
                if best is None or score < best[0]:
                    best = (score, a, b, k, tm)
        return best

    best = search(min_overlap, max_overlap)

    if best is None:
        # No candidate satisfied every target at once (rare, e.g. a squeezed final junction).
        # Fall back to the nominal cut, but still hard-clamp to the 90nt cap - that constraint
        # is never allowed to slip even when Tm/overlap-length can't be perfectly hit here.
        a, b = nominal_start, nominal_end
        if left_start is not None:
            b = min(b, left_start + MAX_OLIGO_LEN)
        if right_end_fixed is not None:
            a = max(a, right_end_fixed - MAX_OLIGO_LEN)
        if b <= a:
            a = b = (a + b) // 2
        return a, b, b - a, calc_tm(sequence[a:b]) if b > a else 0.0
    _, a, b, k, tm = best
    return a, b, k, tm


def build_auto_oligos_pca(sequence, num_oligos, oligo_length, target_tm=TARGET_TM):
    """Like create_balanced_oligos, but nudges each junction to whichever overlap length
    (within MIN_OVERLAP..MAX_OVERLAP, or further if needed) gives a Tm closest to target_tm,
    while keeping every oligo within MIN_OLIGO_LEN..MAX_OLIGO_LEN (90nt cap)."""
    n = len(sequence)
    if num_oligos <= 1 or n <= oligo_length:
        return [(1, n, sequence)]

    nominal = create_balanced_oligos(sequence, num_oligos, oligo_length)
    nominal_bounds = [(s - 1, e) for (s, e, _seq) in nominal]

    starts = [0] * num_oligos
    ends = [0] * num_oligos
    starts[0] = 0
    for j in range(num_oligos - 1):
        nominal_start = nominal_bounds[j + 1][0]
        nominal_end = nominal_bounds[j][1]
        is_last_junction = j == num_oligos - 2
        a, b, _k, _tm = _find_junction_cut(
            sequence, nominal_start, nominal_end, target_tm=target_tm,
            left_start=starts[j], right_end_fixed=n if is_last_junction else None,
            controls_termini=(j % 2 == 0),
        )
        ends[j] = b
        starts[j + 1] = a
    ends[-1] = n

    return [(starts[i] + 1, ends[i], sequence[starts[i]:ends[i]]) for i in range(num_oligos)]


def build_auto_oligos_pca_rtl(sequence, num_oligos, oligo_length, target_tm=TARGET_TM):
    """Same algorithm as build_auto_oligos_pca, but fixes junctions right-to-left instead of
    left-to-right (implemented by running the left-to-right pass on the reversed sequence,
    then mapping the coordinates back). Gives a genuinely different candidate design, since
    whichever end is fixed last is the one most likely to get squeezed."""
    n = len(sequence)
    rev_oligos = build_auto_oligos_pca(sequence[::-1], num_oligos, oligo_length, target_tm)
    oligos = []
    for rev_start, rev_end, rev_seq in rev_oligos:
        orig_start = n - rev_end + 1
        orig_end = n - rev_start + 1
        oligos.append((orig_start, orig_end, rev_seq[::-1]))
    oligos.reverse()
    return oligos


CANDIDATE_STRATEGIES = [
    ("ltr", "왼쪽부터 채우기", None),
    ("rtl", "오른쪽부터 채우기", None),
    ("alt_len", "대체 길이 배치", -2),
]


def build_candidate_oligos(strategy, sequence, num_oligos, oligo_length, target_tm=TARGET_TM):
    key, _label, length_delta = strategy
    length = oligo_length
    if length_delta:
        length = max(MIN_OLIGO_LEN, min(MAX_OLIGO_LEN, oligo_length + length_delta))
    if key == "rtl":
        return build_auto_oligos_pca_rtl(sequence, num_oligos, length, target_tm)
    return build_auto_oligos_pca(sequence, num_oligos, length, target_tm)


UNUSABLE_SCORE = 1_000_000  # disqualifies a candidate from ever ranking first


def score_design(design_summary):
    """Lower is better - aggregate quality score for ranking whole candidate designs.
    A design where any junction's Tm strays more than MAX_TM_DEVIATION from the design's own
    average junction Tm is unusable and is disqualified outright (PCA runs every junction at
    one shared annealing temperature, so a single outlier can break the whole reaction).
    Otherwise, priority order (explicit user preference) is: (1) matching the target Tm and
    being comfortably above TM_COMFORT_MIN (an overlap that only just clears the floor anneals
    less robustly in practice - low fragment concentration, low-GC overlaps, secondary
    structure, and poor transformation efficiency all bite harder at marginal Tm), (2) 3'-end
    quality - dimer/hairpin/mispriming dG and terminal GC clamp, since the last 5-8nt of the 3'
    end matters most for PCR/PCA success and specificity, (3) GC%/repeat-risk/high-T as
    lower-weight tie-breakers, (4) total synthesized nt (cost) as the lowest-priority
    tie-breaker - cheaper is better, but only decides between designs that are otherwise
    equally good; oligo/overlap length itself is unconstrained."""
    s = design_summary
    if not s["usable"]:
        return UNUSABLE_SCORE + s["max_tm_deviation"]
    low_tm_penalty = max(0.0, TM_COMFORT_MIN - s["avg_tm"]) * TM_LOW_PENALTY_WEIGHT if s["avg_tm"] is not None else 0.0
    return (
        s["below_floor"] * 100
        + (s["total_junctions"] - s["sweet"]) * 10
        + low_tm_penalty
        + s.get("dimer_redesign_count", 0) * 80
        + s.get("dimer_caution_count", 0) * 15
        + s.get("terminal_gc_fail_count", 0) * 5
        + len(s["remaining_risk"]) * 4
        + s["high_t_count"] * 3
        + (s["total_junctions"] - s["gc_ok"]) * 2
        + s["total_nt"] * 0.05
    )


def realized_junctions(sequence, oligo_tuples, target_tm=TARGET_TM):
    """Actual overlap length + Tm/GC between each pair of consecutive oligos, from their final boundaries."""
    _, _, floor = tm_bounds(target_tm)
    junctions = []
    for i in range(len(oligo_tuples) - 1):
        end_i = oligo_tuples[i][1]
        start_next = oligo_tuples[i + 1][0] - 1
        if end_i > start_next:
            seg = sequence[start_next:end_i]
            tm = calc_tm(seg)
            junctions.append({
                "overlap_length": end_i - start_next,
                "tm": round(tm, 1),
                "gc_pct": round(calc_gc_fraction(seg) * 100, 0),
                "tm_ok": tm >= floor,
                "risky": _has_homopolymer_run(seg) or _occurs_elsewhere(sequence, seg, start_next),
            })
        else:
            junctions.append({"overlap_length": 0, "tm": None, "gc_pct": None, "tm_ok": None, "risky": False})
    return junctions


def _mark_spans(text, spans):
    """Wrap each (start, end) span of `text` in <mark>, for highlighting risky stretches
    (text is ACGT-only, so this can't inject markup). Overlapping/adjacent spans are merged
    and clamped so the output is always well-formed even if two spans touch."""
    spans = sorted((s for s in spans if s[1] > s[0]), key=lambda s: s[0])
    merged = []
    for s in spans:
        if merged and s[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], s[1]))
        else:
            merged.append(list(s))

    parts = []
    pos = 0
    for start, end in merged:
        start = max(start, pos)
        if start >= end:
            continue
        parts.append(text[pos:start])
        parts.append(f"<mark>{text[start:end]}</mark>")
        pos = end
    parts.append(text[pos:])
    return "".join(parts)


def build_pca_rows(sequence, oligo_tuples, target_tm=TARGET_TM):
    """Alternate strand per oligo (odd index=forward/top strand, even index=reverse/bottom strand)
    and attach the to-next-oligo overlap length + Tm, for a PCA-style annealing assembly.
    The first and last oligo double as the outer PCR primers once PCA is complete."""
    junctions = realized_junctions(sequence, oligo_tuples, target_tm)
    n = len(oligo_tuples)
    rows = []
    for i, (start, end, seq) in enumerate(oligo_tuples):
        is_forward = (i % 2 == 0)
        order_sequence = seq if is_forward else reverse_complement(seq)
        junction = junctions[i] if i < len(junctions) else None
        prev_junction = junctions[i - 1] if i > 0 else None
        t_fraction = order_sequence.count("T") / len(order_sequence) if order_sequence else 0

        # Highlight the risky span(s), if any, in order_sequence coordinates. The overlap with
        # the NEXT oligo sits at the tail of this oligo's top-strand slice; the overlap with the
        # PREVIOUS oligo sits at the head. Reverse-strand oligos read order_sequence backwards,
        # so head/tail swap for them.
        risky_next = bool(junction and junction["risky"])
        risky_prev = bool(prev_junction and prev_junction["risky"])
        oligo_len = len(seq)
        spans = []
        if risky_next:
            k = junction["overlap_length"]
            spans.append((oligo_len - k, oligo_len) if is_forward else (0, k))
        if risky_prev:
            k = prev_junction["overlap_length"]
            spans.append((0, k) if is_forward else (oligo_len - k, oligo_len))
        order_sequence_html = _mark_spans(order_sequence, spans) if spans else order_sequence

        rows.append({
            "index": i + 1,
            "strand": "F" if is_forward else "R",
            "start": start,
            "end": end,
            "length": len(seq),
            "sequence": seq,
            "reverse_complement": reverse_complement(seq),
            "order_sequence": order_sequence,
            "order_sequence_html": order_sequence_html,
            "overlap_to_next": junction["overlap_length"] if junction else None,
            "tm_to_next": junction["tm"] if junction else None,
            "gc_to_next": junction["gc_pct"] if junction else None,
            "tm_ok": junction["tm_ok"] if junction else None,
            "risky_next": risky_next,
            "risky_prev": risky_prev,
            "is_primer": i == 0 or i == n - 1,
            "high_t_warning": t_fraction > HIGH_T_FRACTION,
        })
    return rows


def build_design_summary(sequence, rows, target_tm=TARGET_TM):
    """Short compliance summary against the PCA design rules, re-checked against the final
    (possibly shifted/extended) junctions actually used."""
    if len(rows) < 2:
        return None

    sweet_min, sweet_max, floor = tm_bounds(target_tm)
    within_standard = extended = sweet = floor_only = below_floor = gc_ok = 0
    remaining_risk = []
    for i in range(len(rows) - 1):
        cur, nxt = rows[i], rows[i + 1]
        ov, tm, gc = cur["overlap_to_next"], cur["tm_to_next"], cur["gc_to_next"]
        if ov is None:
            continue
        if ov <= MAX_OVERLAP:
            within_standard += 1
        else:
            extended += 1
        if tm is not None:
            if sweet_min <= tm <= sweet_max:
                sweet += 1
            elif tm >= floor:
                floor_only += 1
            else:
                below_floor += 1
        if gc is not None and GC_MIN * 100 <= gc <= GC_MAX * 100:
            gc_ok += 1

        a, b = nxt["start"] - 1, cur["end"]
        seg = sequence[a:b]
        if _has_homopolymer_run(seg) or _occurs_elsewhere(sequence, seg, a):
            remaining_risk.append(cur["index"])

    total = within_standard + extended
    junction_tms = [r["tm_to_next"] for r in rows[:-1] if r["tm_to_next"] is not None]
    avg_tm = round(sum(junction_tms) / len(junction_tms), 1) if junction_tms else None
    max_tm_deviation = round(max((abs(tm - avg_tm) for tm in junction_tms), default=0.0), 1)
    usable = max_tm_deviation <= MAX_TM_DEVIATION

    overlaps = [r["overlap_to_next"] for r in rows[:-1] if r["overlap_to_next"]]
    overlap_min = min(overlaps) if overlaps else None
    overlap_max = max(overlaps) if overlaps else None
    overlap_avg = round(sum(overlaps) / len(overlaps), 1) if overlaps else None

    return {
        "total_junctions": total,
        "within_standard": within_standard,
        "extended": extended,
        "sweet": sweet,
        "floor_only": floor_only,
        "below_floor": below_floor,
        "gc_ok": gc_ok,
        "remaining_risk": remaining_risk,
        "high_t_count": sum(1 for r in rows if r["high_t_warning"]),
        "target_tm": target_tm,
        "tm_sweet_min": round(sweet_min, 1),
        "tm_sweet_max": round(sweet_max, 1),
        "tm_floor": round(floor, 1),
        "avg_tm": avg_tm,
        "max_tm_deviation": max_tm_deviation,
        "total_nt": sum(r["length"] for r in rows),
        "estimated_price_krw": sum(r["length"] for r in rows) * PRICE_PER_NT_KRW,
        "overlap_min": overlap_min,
        "overlap_max": overlap_max,
        "overlap_avg": overlap_avg,
        "usable": usable,
    }


def pca_rows_to_csv(rows):
    lines = [
        "index,strand,role,start,end,length,order_sequence_5to3,"
        "overlap_to_next_nt,tm_to_next_c,gc_to_next_pct,high_t_warning"
    ]
    for row in rows:
        overlap = row["overlap_to_next"] if row["overlap_to_next"] is not None else ""
        tm = row["tm_to_next"] if row["tm_to_next"] is not None else ""
        gc = row["gc_to_next"] if row["gc_to_next"] is not None else ""
        role = "primer" if row["is_primer"] else ""
        lines.append(
            f"{row['index']},{row['strand']},{role},{row['start']},{row['end']},{row['length']},"
            f"{row['order_sequence']},{overlap},{tm},{gc},{row['high_t_warning']}"
        )
    return "\n".join(lines) + "\n"


def pca_rows_to_fasta(rows):
    lines = []
    for row in rows:
        lines.append(f">oligo_{row['index']}_{row['strand']}_{row['start']}-{row['end']}")
        lines.append(row["order_sequence"])
    return "\n".join(lines) + "\n"


def build_diagram_svg(rows, total_length, width=860, pad=20):
    """A two-lane schematic: forward oligos on the top lane, reverse oligos on the bottom
    lane, both positioned along the true sequence coordinate, with overlap + Tm labels
    connecting consecutive oligos."""
    if not rows or total_length <= 0:
        return ""

    usable = width - 2 * pad
    scale = usable / total_length
    row_h = 22
    top_y = 28
    bottom_y = 96
    height = bottom_y + row_h + 34

    def x_of(pos0):
        return pad + pos0 * scale

    colors = {"F": "#2563eb", "R": "#dc2626"}
    parts = [
        f'<svg viewBox="0 0 {width} {height}" width="100%" xmlns="http://www.w3.org/2000/svg" '
        f'font-family="ui-monospace, Menlo, Consolas, monospace" font-size="10">',
        f'<line x1="{pad}" y1="{top_y + row_h / 2:.1f}" x2="{width - pad}" y2="{top_y + row_h / 2:.1f}" '
        f'stroke="#9aa1ac" stroke-width="1" stroke-dasharray="2,3" />',
        f'<line x1="{pad}" y1="{bottom_y + row_h / 2:.1f}" x2="{width - pad}" y2="{bottom_y + row_h / 2:.1f}" '
        f'stroke="#9aa1ac" stroke-width="1" stroke-dasharray="2,3" />',
        f'<text x="{pad}" y="{top_y - 8}" fill="#6b7280">Forward (top strand, 5\'&#8594;3\')</text>',
        f'<text x="{pad}" y="{bottom_y + row_h + 14}" fill="#6b7280">Reverse (bottom strand, 5\'&#8594;3\')</text>',
        f'<text x="{width - pad}" y="{bottom_y + row_h + 14}" text-anchor="end" fill="#f59e0b">'
        f'&#9633; amber outline = also used as PCR primer</text>',
    ]

    for row in rows:
        y = top_y if row["strand"] == "F" else bottom_y
        x1 = x_of(row["start"] - 1)
        w = max((row["end"] - row["start"] + 1) * scale, 2)
        color = colors[row["strand"]]
        parts.append(f'<rect x="{x1:.1f}" y="{y}" width="{w:.1f}" height="{row_h}" rx="3" fill="{color}" opacity="0.88" />')
        if row["is_primer"]:
            parts.append(
                f'<rect x="{x1:.1f}" y="{y}" width="{w:.1f}" height="{row_h}" rx="3" '
                f'fill="none" stroke="#f59e0b" stroke-width="2" />'
            )
        arrow = "&#9654;" if row["strand"] == "F" else "&#9664;"
        label = f'{arrow} {row["strand"]}{row["index"]}' if row["strand"] == "R" else f'{row["strand"]}{row["index"]} {arrow}'
        parts.append(
            f'<text x="{x1 + w / 2:.1f}" y="{y + row_h / 2 + 3.5:.1f}" text-anchor="middle" '
            f'fill="white" font-weight="700">{label}</text>'
        )

    for i in range(len(rows) - 1):
        cur, nxt = rows[i], rows[i + 1]
        if not cur["overlap_to_next"]:
            continue
        mid_x = x_of(((nxt["start"] - 1) + cur["end"]) / 2)
        line_color = "#dc2626" if cur["tm_ok"] is False else "#9aa1ac"
        parts.append(
            f'<line x1="{mid_x:.1f}" y1="{top_y + row_h}" x2="{mid_x:.1f}" y2="{bottom_y}" '
            f'stroke="{line_color}" stroke-width="1.5" />'
        )
        tm_txt = f'{cur["overlap_to_next"]}nt / {cur["tm_to_next"]}&#176;C' if cur["tm_to_next"] is not None else f'{cur["overlap_to_next"]}nt'
        parts.append(
            f'<text x="{mid_x:.1f}" y="{(top_y + row_h + bottom_y) / 2 + 3.5:.1f}" text-anchor="middle" '
            f'fill="{line_color if cur["tm_ok"] is False else "#6b7280"}">{tm_txt}</text>'
        )

    parts.append("</svg>")
    return "".join(parts)


def tm_scan_values():
    values = []
    t = TM_SCAN_MIN
    while t <= TM_SCAN_MAX + 1e-9:
        values.append(round(t, 1))
        t += TM_SCAN_STEP
    return values


def scan_auto_candidates(sequence, candidate_key=None):
    """Shared by Oligo Assembly Designer and Silent Library Designer's auto mode: plan the
    oligo count/length from `sequence`'s length, scan the whole practical Tm range (see
    TM_SCAN_MIN/MAX/STEP) x 3 build strategies, run 3'-end dimer/hairpin/GC-clamp QC on each,
    and rank every combination (see score_design). Both tools apply this to the same real,
    fixed DNA sequence - Silent Library Designer's IUPAC-degenerate positions have no defined
    Tm, but the oligo junctions themselves must still anneal reliably, so boundaries are
    always chosen against the real sequence.

    Returns (candidates, primary, computed_oligo_length, error). `candidates` is capped at
    MAX_CANDIDATES_SHOWN and sorted best-first; `primary` is the requested `candidate_key` or
    the best one. On failure, candidates/primary/computed_oligo_length are None and `error`
    is set - either planning failed outright, or (when candidate_key wasn't explicitly given)
    even the best candidate is unusable.
    """
    plan = plan_oligo_split(len(sequence), NOMINAL_OVERLAP_FOR_PLANNING)
    if plan is None:
        error = (
            f"자동 설계 실패: 서열 길이({len(sequence)}nt)로는 올리고 길이를 "
            f"{MIN_OLIGO_LEN}~{MAX_OLIGO_LEN}nt 범위 안에서 짝수 개로 나눌 수 없습니다. "
            "수동 모드를 사용해보세요."
        )
        return None, None, None, error
    num_oligos, computed_oligo_length = plan

    per_tm_best = []
    for trial_tm in tm_scan_values():
        best_for_tm = None
        for strategy in CANDIDATE_STRATEGIES:
            key, label, _length_delta = strategy
            oligo_tuples = build_candidate_oligos(strategy, sequence, num_oligos, computed_oligo_length, trial_tm)
            cand_rows = build_pca_rows(sequence, oligo_tuples, trial_tm)
            cand_summary = build_design_summary(sequence, cand_rows, trial_tm)
            score = score_design(cand_summary) if cand_summary else float("inf")
            if best_for_tm is None or score < best_for_tm["score"]:
                best_for_tm = {
                    "key": f"{key}@{trial_tm}",
                    "label": f"{label} (Tm {trial_tm}°C)",
                    "target_tm": trial_tm,
                    "rows": cand_rows,
                    "total_nt": sum(row["length"] for row in cand_rows),
                    "estimated_price_krw": sum(row["length"] for row in cand_rows) * PRICE_PER_NT_KRW,
                    "design_summary": cand_summary,
                    "score": score,
                }
        per_tm_best.append(best_for_tm)

    for c in per_tm_best:
        add_primer_qc(c["rows"], sequence)
        s = c["design_summary"]
        if s is not None:
            s["dimer_redesign_count"] = sum(1 for r in c["rows"] if r["dimer_flag"] == "redesign")
            s["dimer_caution_count"] = sum(1 for r in c["rows"] if r["dimer_flag"] == "caution")
            s["terminal_gc_fail_count"] = sum(1 for r in c["rows"] if not r["terminal_gc_ok"])
            c["score"] = score_design(s)

    all_sorted = sorted(per_tm_best, key=lambda c: c["score"])
    best = all_sorted[0]
    if candidate_key is None and best["design_summary"] and not best["design_summary"]["usable"]:
        dev = best["design_summary"]["max_tm_deviation"]
        error = (
            f"사용 가능한 설계를 찾지 못했습니다: {TM_SCAN_MIN}~{TM_SCAN_MAX}°C 전체를 탐색했지만, "
            f"최선의 결과도 접합부 간 Tm 편차가 최대 {dev}°C로 허용 기준"
            f"(평균 대비 {MAX_TM_DEVIATION}°C)을 초과합니다."
        )
        return None, None, None, error

    candidates = all_sorted[:MAX_CANDIDATES_SHOWN]
    for i, c in enumerate(candidates):
        c["rank"] = i + 1
    if candidate_key and not any(c["key"] == candidate_key for c in candidates):
        selected = next((c for c in all_sorted if c["key"] == candidate_key), None)
        if selected:
            candidates = candidates[:-1] + [selected]

    primary = next((c for c in candidates if c["key"] == candidate_key), best)
    return candidates, primary, computed_oligo_length, None
