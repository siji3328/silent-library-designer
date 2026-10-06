# -*- coding: utf-8 -*-
"""Silent Library Designer - synonymous (silent) mutation oligo library design."""

from collections import Counter

from .common import (
    TARGET_TM,
    _mark_spans,
    add_primer_qc,
    build_design_summary,
    build_diagram_svg,
    build_pca_rows,
    create_oligos,
    reverse_complement,
    scan_auto_candidates,
)

STANDARD_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

IUPAC_CODON_TABLE = {
    "F": "TTY", "L": "CTN", "I": "ATH", "M": "ATG", "V": "GTN",
    "S": "TCN", "P": "CCN", "T": "ACN", "A": "GCN", "Y": "TAY",
    "H": "CAY", "Q": "CAR", "N": "AAY", "K": "AAR", "D": "GAY",
    "E": "GAR", "C": "TGY", "W": "TGG", "R": "CGN", "G": "GGN",
    "*": "TAR",
}

IUPAC_DEGENERATE_CODES = "NRYMKSWHBVD"

VALID_AA_LETTERS = set("ACDEFGHIKLMNPQRSTVWY*")

DEFAULT_DNA_SEQUENCE = (
    "ATGAGCAAAGGTGAAGAACTGTTTACCGGCGTTGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAACGGTCACAAATTCAGC"
    "GTGCGTGGTGAAGGTGAAGGCGATGCCACGATTGGCAAACTGACGCTGAAATTTATCTGCACCACCGGCAAACTGCCGGTGCCG"
    "TGGCCGACGCTGGTGACCACCCTGACCTATGGCGTTCAGTGTTTTAGTCGCTATCCGGATCACATGAAACGTCACGATTTCTTT"
    "AAATCTGCAATGCCGGAAGGCTATGTGCAGGAACGTACGATTAGCTTTAAAGATGATGGCAAATATAAAACGCGCGCCGTTGTG"
    "AAATTTGAAGGCGATACCCTGGTGAACCGCATTGAACTGAAAGGCACGGATTTTAAAGAAGATGGCAATATCCTGGGCCATAAA"
    "CTGGAATACAACTTTAATAGCCATAATGTTTATATTACGGCGGATAAACAGAAAAATGGCATCAAAGCGAATTTTACCGTTCGC"
    "CATAACGTTGAAGATGGCAGTGTGCAGCTGGCAGATCATTATCAGCAGAATACCCCGATTGGTGATGGTCCGGTGCTGCTGCCG"
    "GATAATCATTATCTGAGCACGCAGACCGTTCTGTCTAAAGATCCGAACGAAAAACGGGACCACATGGTTCTGCACGAATATGTG"
    "AATGCGGCAGGTATTACGTAA"
)

DEFAULT_AA_SEQUENCE = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATIGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFF"
    "KSAMPEGYVQERTISFKDDGKYKTRAVVKFEGDTLVNRIELKGTDFKEDGNILGHKLEYNFNSHNVYITADKQKNGIKANFTVR"
    "HNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQTVLSKDPNEKRDHMVLHEYVNAAGIT*"
)


def get_iupac_codon(codon, amino_acid):
    """Pick the degenerate (IUPAC) codon that stays synonymous with the observed codon."""
    if amino_acid == "L":
        return "CTN" if codon.startswith("CT") else "TTR"
    elif amino_acid == "S":
        return "TCN" if codon.startswith("TC") else "AGY"
    elif amino_acid == "R":
        return "CGN" if codon.startswith("CG") else "AGR"
    elif amino_acid == "*":
        return "TAR" if codon.startswith("TA") else "TGA"
    else:
        return IUPAC_CODON_TABLE.get(amino_acid, codon)


def translate(dna_sequence):
    protein = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i + 3]
        protein.append(STANDARD_CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def build_iupac_sequence(dna_sequence, amino_acid_sequence):
    parts = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        amino_acid = amino_acid_sequence[i // 3]
        parts.append(get_iupac_codon(codon, amino_acid))
    return "".join(parts)


def fix_overlap_codons(dna_sequence, iupac_sequence, oligo_tuples):
    """The region two neighboring oligos share has to physically anneal, so it can't be
    independently degenerate on each oligo - two separately-synthesized degenerate pools
    would not reliably end up complementary there. Patch every junction's overlap span
    (extended to whole-codon boundaries, so no codon is left half-fixed/half-degenerate) back
    to the original, non-degenerate DNA codons.

    Returns (patched_sequence, fixed_ranges) where fixed_ranges is a list of (start0, end0)
    codon-aligned spans that lost degeneracy because of this.
    """
    chars = list(iupac_sequence)
    fixed_ranges = []
    for i in range(len(oligo_tuples) - 1):
        end_i0 = oligo_tuples[i][1]
        next_start0 = oligo_tuples[i + 1][0] - 1
        if end_i0 <= next_start0:
            continue
        ext_start = (next_start0 // 3) * 3
        ext_end = min(-(-end_i0 // 3) * 3, len(chars))
        if any(chars[p] != dna_sequence[p] for p in range(ext_start, ext_end)):
            fixed_ranges.append((ext_start, ext_end))
        for pos in range(ext_start, ext_end):
            chars[pos] = dna_sequence[pos]
    return "".join(chars), fixed_ranges


def _highlight_fixed(sequence, fixed_ranges):
    if not fixed_ranges:
        return sequence
    parts = []
    pos = 0
    for start, end in sorted(fixed_ranges):
        if start > pos:
            parts.append(sequence[pos:start])
        parts.append(f'<span class="fixed-codon">{sequence[start:end]}</span>')
        pos = max(pos, end)
    parts.append(sequence[pos:])
    return "".join(parts)


def apply_library_content(rows, dna_sequence, iupac_sequence):
    """Junction boundaries are chosen against the real, fixed DNA sequence (see
    scan_auto_candidates / common.build_pca_rows) so Tm/GC/dimer checks stay meaningful. This
    re-derives each oligo's actual sequence for the library: degenerate (IUPAC) in its own
    unique body, but fixed back to the original codons across every junction overlap (see
    fix_overlap_codons) so neighboring oligos still reliably anneal there.

    Mutates `rows` in place (sequence/reverse_complement/order_sequence/order_sequence_html)
    and returns (final_sequence, fixed_ranges).
    """
    oligo_tuples = [(r["start"], r["end"], dna_sequence[r["start"] - 1:r["end"]]) for r in rows]
    final_sequence, fixed_ranges = fix_overlap_codons(dna_sequence, iupac_sequence, oligo_tuples)

    for i, row in enumerate(rows):
        seq = final_sequence[row["start"] - 1:row["end"]]
        is_forward = row["strand"] == "F"
        order_sequence = seq if is_forward else reverse_complement(seq)
        oligo_len = len(seq)

        spans = []
        if row.get("risky_next"):
            k = row["overlap_to_next"]
            spans.append((oligo_len - k, oligo_len) if is_forward else (0, k))
        if row.get("risky_prev"):
            prev_overlap = rows[i - 1]["overlap_to_next"] if i > 0 else None
            if prev_overlap:
                spans.append((0, prev_overlap) if is_forward else (oligo_len - prev_overlap, oligo_len))

        row["sequence"] = seq
        row["reverse_complement"] = reverse_complement(seq)
        row["order_sequence"] = order_sequence
        row["order_sequence_html"] = _mark_spans(order_sequence, spans) if spans else order_sequence

    return final_sequence, fixed_ranges


def run(dna_sequence, amino_acid_sequence, oligo_length, overlap, auto_split, candidate_key=None):
    """Validate inputs and compute the silent-mutation library + oligo split.

    Junction boundaries (and, in auto mode, the whole Tm scan/candidate ranking) are computed
    against the real, fixed `dna_sequence` - identical to Oligo Assembly Designer - since the
    physical annealing at each overlap depends on the real bases there, not the degenerate
    library encoding. The IUPAC degenerate codons are layered back in afterward, everywhere
    except the overlaps themselves (see apply_library_content).

    Returns (result_dict_or_None, error_or_None, warning_or_None).
    """
    if not dna_sequence:
        return None, "DNA 서열을 입력해주세요.", None
    if any(b not in "ATGC" for b in dna_sequence):
        return None, "DNA 서열에는 A, T, G, C만 사용할 수 있습니다.", None
    if not auto_split and oligo_length <= overlap:
        return None, "올리고 길이는 오버랩 길이보다 커야 합니다.", None
    if not amino_acid_sequence:
        return None, "아미노산(AA) 서열을 입력해주세요.", None
    if len(dna_sequence) % 3 != 0:
        return None, f"DNA 서열 길이({len(dna_sequence)})가 3의 배수가 아닙니다.", None
    if any(a not in VALID_AA_LETTERS for a in amino_acid_sequence):
        return None, "AA 서열에는 표준 아미노산 1글자 코드(및 정지코돈 '*')만 사용할 수 있습니다.", None
    if len(dna_sequence) != len(amino_acid_sequence) * 3:
        return None, (
            f"DNA 서열 길이({len(dna_sequence)})와 AA 서열 길이({len(amino_acid_sequence)}"
            f" x 3 = {len(amino_acid_sequence) * 3})가 일치하지 않습니다."
        ), None

    warning = None
    translated = translate(dna_sequence)
    if translated != amino_acid_sequence:
        mismatches = [
            i + 1
            for i, (a, b) in enumerate(zip(translated, amino_acid_sequence))
            if a != b
        ]
        warning = (
            "입력하신 DNA 서열을 번역한 결과가 AA 서열과 다릅니다 "
            f"(불일치 위치: {', '.join(map(str, mismatches[:10]))}"
            f"{' 외 다수' if len(mismatches) > 10 else ''}). "
            "AA 서열을 기준으로 계산을 진행합니다."
        )

    iupac_sequence = build_iupac_sequence(dna_sequence, amino_acid_sequence)

    computed_oligo_length = None
    candidates = None

    if auto_split:
        candidates, primary, computed_oligo_length, error = scan_auto_candidates(dna_sequence, candidate_key)
        if error:
            return None, error, warning
        rows = primary["rows"]
        design_summary = primary["design_summary"]
    else:
        oligo_tuples = create_oligos(dna_sequence, oligo_length, overlap)
        if len(oligo_tuples) % 2 != 0:
            warning = (
                (warning + " ") if warning else ""
            ) + (
                f"올리고 개수가 홀수({len(oligo_tuples)}개)입니다. Forward/Reverse가 깔끔하게 "
                "번갈아 이어지려면 짝수 개를 권장합니다 — 올리고 길이/오버랩을 조정하거나 자동 최적화를 사용해보세요."
            )
        rows = build_pca_rows(dna_sequence, oligo_tuples, TARGET_TM)
        design_summary = build_design_summary(dna_sequence, rows, TARGET_TM)
        add_primer_qc(rows, dna_sequence)
        if design_summary is not None:
            design_summary["dimer_redesign_count"] = sum(1 for r in rows if r["dimer_flag"] == "redesign")
            design_summary["dimer_caution_count"] = sum(1 for r in rows if r["dimer_flag"] == "caution")
            design_summary["terminal_gc_fail_count"] = sum(1 for r in rows if not r["terminal_gc_ok"])

    # Junction boundaries above were chosen against the real, fixed DNA sequence so Tm/GC/
    # dimer analysis stays meaningful. Now layer the IUPAC degenerate codons back in for the
    # actual library - everywhere except each junction's overlap, which must stay fixed so
    # neighboring oligos still reliably anneal there (see apply_library_content).
    final_sequence, fixed_ranges = apply_library_content(rows, dna_sequence, iupac_sequence)
    diagram_svg = build_diagram_svg(rows, len(dna_sequence))

    counter = Counter(c for c in final_sequence if c in IUPAC_DEGENERATE_CODES)
    degenerate_counts = [(code, counter[code]) for code in IUPAC_DEGENERATE_CODES]
    fixed_codon_count = sum((end - start) // 3 for start, end in fixed_ranges)

    result = {
        "amino_acid_sequence": amino_acid_sequence,
        "iupac_sequence": final_sequence,
        "iupac_sequence_html": _highlight_fixed(final_sequence, fixed_ranges),
        "fixed_codon_count": fixed_codon_count,
        "dna_length": len(dna_sequence),
        "degenerate_counts": degenerate_counts,
        "oligos": rows,
        "auto_split": auto_split,
        "computed_oligo_length": computed_oligo_length,
        "diagram_svg": diagram_svg,
        "design_summary": design_summary,
        "candidates": candidates,
        "primary_key": primary["key"] if candidates else None,
    }
    return result, None, warning
