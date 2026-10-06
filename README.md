---
title: Silent Library Designer
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 7860
pinned: false
license: mit
short_description: Synonymous DNA library oligomer designer with assembly-aware overlaps
---

# Silent Library Designer

**Synonymous DNA Library Oligomer Designer**

Design assembly-ready oligomers for synonymous (silent) DNA libraries from a wild-type
coding sequence. Paste a CDS and its amino-acid sequence; the tool chooses oligomer
junctions, applies degenerate codons, and returns an order-ready oligomer table.

> 야생형 코딩 서열(CDS)을 넣으면, 아미노산 서열은 그대로 두고 동의 코돈만 바꾼 라이브러리를
> 만든 뒤 PCA 조립이 가능하도록 올리고머로 나눠 주는 웹 도구입니다.

| | |
|---|---|
| **Input** | WT coding sequence + its amino-acid sequence (length ÷ 3 must match) |
| **Output** | IUPAC degenerate library sequence, split into PCA-ready oligomers, with per-junction Tm/GC/dimer QC and CSV/FASTA export |

## How the design works

The order of operations is the point, and it is not the obvious one.

1. **Junctions are chosen on the real WT bases first.** Nearest-neighbor Tm, GC content,
   3′ terminal clamp, and dimer/hairpin ΔG are all computed on actual bases. In auto mode
   the tool sweeps 45–65 °C across three assembly strategies, scores every candidate
   design, and keeps the best. A design where any junction's Tm strays too far from that
   design's own average is disqualified outright, since PCA runs every junction at one
   shared annealing temperature.
2. **Degenerate codons are layered in afterwards.** Each codon becomes the IUPAC code that
   stays synonymous given the observed codon — the branch matters, so Leu is `CTN` or
   `TTR` depending on its first two bases, and likewise for Ser, Arg and stop.
3. **Every junction overlap is reverted to the original codons.** Two independently
   synthesized degenerate pools would not reliably be complementary across an overlap, so
   annealing there could not be counted on and the resulting library composition would be
   unpredictable. The overlaps stay fixed; the tool reports how many codons lost
   degeneracy as a result.

Because the overlaps are fixed, the Tm computed in step 1 is the Tm of the molecules that
are actually ordered. The cost is that those codons carry no library diversity, which the
result page states explicitly.

## Running locally

Only Flask is required — the thermodynamics are implemented in `tools/common.py` with no
external scientific dependencies.

```bash
pip install -r requirements.txt
python app.py
```

Then open <http://127.0.0.1:7860> — the same port the container serves, so local and
deployed runs agree. Set `PORT` to override it.

## Layout

```
app.py                      routing and form handling
tools/common.py             Tm / GC / dimer scoring, junction search, oligomer rows
tools/silent_library.py     degenerate codon table, overlap fixing
templates/index.html        the single page
static/style.css            stylesheet
static/fonts/               bundled DejaVu Sans / Sans Mono
Dockerfile                  Hugging Face Spaces (Docker SDK, port 7860)
```

## Notes

For research use only. Check synthesized sequences for restriction sites and unwanted
repeats before ordering.

Bundled DejaVu fonts are distributed under the [DejaVu Fonts License](https://dejavu-fonts.github.io/License.html).
