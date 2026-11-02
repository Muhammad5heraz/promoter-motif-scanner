# promoter-motif-scanner

A small, dependency-free Python tool for scanning promoter sequences for
regulatory motifs. It reports every hit on both strands, summarises counts per
sequence, and flags sequences where two motifs occur within a chosen distance
of each other (for example a GA-repeat and a C-box / E-box element).

Built from the kind of promoter analysis used in plant transcription-factor
gene-family studies, written here as a general, reusable script.

## Features
- IUPAC-coded motifs (`R`, `Y`, `N`, ...), so degenerate sites work
- Both strands, with palindromic motifs counted once
- Overlapping matches reported (e.g. `GAGAGA` inside `GAGAGAGA`)
- 1-based coordinates in the forward-strand frame
- Optional per-sequence summary and two-motif co-occurrence flag
- Standard library only (Python 3.8+)

## Usage
```bash
python motif_scan.py example/promoters.fa \
  -m GA_repeat=GAGAGA -m Ebox=CACGTG \
  -w 50 -s summary.tsv > hits.tsv
```

| Option | Meaning |
|---|---|
| `fasta` | Input promoter FASTA |
| `-m NAME=SEQ` | Motif (IUPAC); repeat for several. `SEQ` alone also works |
| `-o FILE` | Hits table (default: stdout) |
| `-s FILE` | Per-sequence summary table |
| `-w N` | Co-occurrence window in bp for the first two motifs (default 100) |

## Output
`hits.tsv`: `seq_id  motif  strand  start  end  match`

`summary.tsv`: `seq_id  length  <motif>_hits ...  cooccur_within_<N>bp`

A count of sequences with at least one hit per motif, and with the two
motifs co-occurring, is printed to stderr.

## Example
Running the command above on `example/promoters.fa` (three toy promoters)
finds GA-repeat hits on both strands and reports two sequences where the
GA-repeat and the E-box fall within 50 bp.

## Tests
```bash
python tests/test_scanner.py
```

## Notes
This is a scanning utility. A hit is a sequence match, not evidence of
binding or regulation; enrichment should be assessed against a suitable
background (for example a genome-wide scan) before drawing conclusions.

## Author
Muhammad Sheraz · [ORCID 0009-0007-3877-2232](https://orcid.org/0009-0007-3877-2232)

## License
MIT
