#!/usr/bin/env python3
"""promoter-motif-scanner

Scan promoter sequences (FASTA) for one or more IUPAC-coded motifs on both
strands, report every hit with 1-based coordinates, and summarise per-sequence
counts and motif co-occurrence within a distance window.

Pure standard library; no dependencies.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}
COMP = {
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M",
    "M": "K", "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
}


def iupac_to_regex(motif):
    """Compile an IUPAC motif into a regex (overlapping matches via lookahead)."""
    try:
        body = "".join(IUPAC[c] for c in motif.upper())
    except KeyError as bad:
        raise ValueError(f"Invalid IUPAC code {bad} in motif '{motif}'")
    return re.compile(f"(?=({body}))")


def revcomp_motif(motif):
    return "".join(COMP[c] for c in reversed(motif.upper()))


def read_fasta(handle):
    """Yield (id, sequence) from a FASTA file handle."""
    name, chunks = None, []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def scan_sequence(seq, motif_name, motif):
    """Return hits as (strand, start, end, matched_text), 1-based inclusive."""
    hits = []
    length = len(motif)
    fwd = iupac_to_regex(motif)
    for m in fwd.finditer(seq):
        s = m.start()
        hits.append(("+", s + 1, s + length, m.group(1)))
    rc = revcomp_motif(motif)
    if rc != motif.upper():  # palindromes would be double counted
        rev = iupac_to_regex(rc)
        for m in rev.finditer(seq):
            s = m.start()
            hits.append(("-", s + 1, s + length, m.group(1)))
    return sorted(hits, key=lambda h: (h[1], h[0]))


def cooccurs(hits_a, hits_b, window):
    """True if any hit of A lies within `window` bp of any hit of B."""
    for _, sa, ea, _ in hits_a:
        for _, sb, eb, _ in hits_b:
            gap = max(sb - ea, sa - eb, 0)  # 0 if overlapping
            if gap <= window:
                return True
    return False


def parse_motifs(items):
    motifs = []
    for item in items:
        if "=" in item:
            name, seq = item.split("=", 1)
        else:
            name, seq = item, item
        motifs.append((name, seq))
    return motifs


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Scan promoter FASTA for IUPAC motifs on both strands.")
    ap.add_argument("fasta", help="Input promoter FASTA")
    ap.add_argument("-m", "--motif", action="append", required=True,
                    help="Motif as NAME=SEQUENCE or SEQUENCE (IUPAC); repeatable")
    ap.add_argument("-o", "--hits", default="-",
                    help="Hits TSV output (default: stdout)")
    ap.add_argument("-s", "--summary", help="Per-sequence summary TSV")
    ap.add_argument("-w", "--window", type=int, default=100,
                    help="Co-occurrence distance in bp for the first two "
                         "motifs (default: 100)")
    args = ap.parse_args(argv)

    motifs = parse_motifs(args.motif)
    out = sys.stdout if args.hits == "-" else open(args.hits, "w", newline="")
    hit_writer = csv.writer(out, delimiter="\t")
    hit_writer.writerow(["seq_id", "motif", "strand", "start", "end", "match"])

    summary_rows, totals = [], defaultdict(int)
    n_seq = n_co = 0
    with open(args.fasta) as fh:
        for seq_id, seq in read_fasta(fh):
            n_seq += 1
            per_motif = {}
            for name, motif in motifs:
                hits = scan_sequence(seq, name, motif)
                per_motif[name] = hits
                for strand, s, e, text in hits:
                    hit_writer.writerow([seq_id, name, strand, s, e, text])
                if hits:
                    totals[name] += 1
            co = ""
            if len(motifs) >= 2:
                a, b = motifs[0][0], motifs[1][0]
                co = int(cooccurs(per_motif[a], per_motif[b], args.window))
                n_co += co
            summary_rows.append(
                [seq_id, len(seq)] + [len(per_motif[n]) for n, _ in motifs] + [co])
    if out is not sys.stdout:
        out.close()

    if args.summary:
        with open(args.summary, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            head = ["seq_id", "length"] + [f"{n}_hits" for n, _ in motifs]
            w.writerow(head + [f"cooccur_within_{args.window}bp"])
            w.writerows(summary_rows)

    print(f"Sequences scanned: {n_seq}", file=sys.stderr)
    for name, _ in motifs:
        print(f"  {name}: {totals[name]} sequences with >=1 hit", file=sys.stderr)
    if len(motifs) >= 2:
        print(f"  {motifs[0][0]} + {motifs[1][0]} within {args.window} bp: "
              f"{n_co} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
