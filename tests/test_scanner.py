import io
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from motif_scan import (cooccurs, read_fasta, revcomp_motif,  # noqa: E402
                        scan_sequence)


def test_revcomp_iupac():
    assert revcomp_motif("GAGAGA") == "TCTCTC"
    assert revcomp_motif("CANNTG") == "CANNTG"  # E-box style palindrome


def test_both_strands_and_overlap():
    hits = scan_sequence("GAGAGAGA", "ga", "GAGA")
    plus = [h for h in hits if h[0] == "+"]
    assert [h[1] for h in plus] == [1, 3, 5]  # overlapping matches kept
    minus = scan_sequence("TCTCTC", "ga", "GAGA")
    assert any(h[0] == "-" for h in minus)


def test_palindrome_not_double_counted():
    hits = scan_sequence("AACGTTT", "cg", "ACGT")
    assert len(hits) == 1


def test_cooccurrence_window():
    a = [("+", 10, 15, "x")]
    b = [("+", 60, 65, "y")]
    assert cooccurs(a, b, 50)
    assert not cooccurs(a, b, 40)


def test_fasta_reader():
    seqs = list(read_fasta(io.StringIO(">s1 desc\nAC\nGT\n>s2\nTT\n")))
    assert seqs == [("s1", "ACGT"), ("s2", "TT")]


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
