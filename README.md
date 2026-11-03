# promoter-motif-scanner

Scan promoter sequences for IUPAC-coded regulatory motifs on both strands and
flag sequences where two motifs occur close together.

The tool reports every hit with 1-based coordinates, summarises hit counts per
sequence, and marks sequences in which the first two motifs fall within a
chosen distance of each other (for example a GA-repeat and a G-box). It is a
general, reusable script of the kind used for promoter analysis in plant
transcription-factor gene-family studies.

## Features
- IUPAC-coded motifs (`R`, `Y`, `N`, ...), so degenerate sites work
- Both strands scanned; palindromic motifs (e.g. `CACGTG`) counted once
- Overlapping matches reported (e.g. `GAGAGA` three times inside `GAGAGAGAGA`)
- 1-based coordinates in the forward-strand frame
- Optional per-sequence summary and two-motif co-occurrence flag
- Standard library only

## Requirements and installation
- Python 3.8 or later
- No third-party packages

```bash
git clone https://github.com/Muhammad5heraz/promoter-motif-scanner.git
cd promoter-motif-scanner
python motif_scan.py --help
```

## Usage
```bash
python motif_scan.py example/promoters.fa \
  -m GA_repeat=GAGAGA -m G_box=CACGTG \
  -w 50 -s summary.tsv > hits.tsv
```

| Option | Meaning |
|---|---|
| `fasta` | Input promoter FASTA |
| `-m NAME=SEQ` | Motif in IUPAC code; repeat for several motifs. `SEQ` alone also works and is used as the name |
| `-o FILE` | Hits table (default: stdout) |
| `-s FILE` | Per-sequence summary table |
| `-w N` | Co-occurrence window in bp, applied to the first two motifs (default 100) |

## Output
**Hits table**: one row per hit.

`seq_id  motif  strand  start  end  match`

On the `-` strand, `match` is the forward-strand text of the site (the reverse
complement of the motif), and `start`/`end` are forward-strand coordinates.

**Summary table** (`-s`): one row per sequence.

`seq_id  length  <motif>_hits ...  cooccur_within_<N>bp`

The co-occurrence column (`1` = yes, `0` = no) is written only when two or more
motifs are given.

**stderr**: the number of sequences scanned, the number with at least one hit
per motif, and the number with the first two motifs co-occurring.

## Example
Running the usage command above on `example/promoters.fa` (three short toy
promoters) produces:

`hits.tsv`
```
seq_id          motif      strand  start  end  match
geneA_promoter  GA_repeat  +       4      9    GAGAGA
geneA_promoter  GA_repeat  +       6      11   GAGAGA
geneA_promoter  GA_repeat  +       8      13   GAGAGA
geneA_promoter  GA_repeat  +       40     45   GAGAGA
geneA_promoter  G_box      +       22     27   CACGTG
geneC_promoter  GA_repeat  -       3      8    TCTCTC
geneC_promoter  GA_repeat  -       5      10   TCTCTC
geneC_promoter  GA_repeat  -       7      12   TCTCTC
geneC_promoter  GA_repeat  +       40     45   GAGAGA
geneC_promoter  GA_repeat  +       42     47   GAGAGA
geneC_promoter  GA_repeat  +       44     49   GAGAGA
geneC_promoter  G_box      +       19     24   CACGTG
```

`summary.tsv`
```
seq_id          length  GA_repeat_hits  G_box_hits  cooccur_within_50bp
geneA_promoter  47      4               1           1
geneB_promoter  38      0               0           0
geneC_promoter  49      6               1           1
```

stderr
```
Sequences scanned: 3
  GA_repeat: 2 sequences with >=1 hit
  G_box: 2 sequences with >=1 hit
  GA_repeat + G_box within 50 bp: 2 sequences
```

(Columns are tab-separated in the real files; they are aligned here for
reading.)

## Tests
```bash
python tests/test_scanner.py
```
Six tests cover reverse complementation of IUPAC codes, both-strand and
overlapping matches, palindrome handling, the co-occurrence window, FASTA
parsing and the single-motif summary layout.

## Limitations
- A hit is an exact sequence match, not evidence of protein binding or
  regulation. Assess enrichment against a suitable background (for example a
  genome-wide promoter scan) before drawing conclusions.
- Co-occurrence is computed for the first two motifs only.
- No position weight matrices or mismatch tolerance; use dedicated tools
  (for example FIMO in the MEME Suite) for those.
- Motif names should be unique; a repeated name overwrites the earlier motif
  in the summary.

## How to cite
If you use this tool, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button).

## Author
Muhammad Sheraz · [ORCID 0009-0007-3877-2232](https://orcid.org/0009-0007-3877-2232)

## License
MIT. See [`LICENSE`](LICENSE).
