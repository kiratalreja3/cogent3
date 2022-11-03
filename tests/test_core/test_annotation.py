#!/usr/bin/env python

import unittest

from cogent3 import DNA, make_aligned_seqs
from cogent3.core.annotation import Feature, _Feature
from cogent3.core.location import Map, Span, as_map
from cogent3.core.sequence import DnaSequence, RnaSequence


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


def makeSampleSequence(with_gaps=False):
    raw_seq = "AACCCAAAATTTTTTGGGGGGGGGGCCCC"
    cds = (15, 25)
    utr = (12, 15)
    if with_gaps:
        raw_seq = raw_seq[:5] + "-----" + raw_seq[10:-2] + "--"
    seq = DNA.make_seq(raw_seq)
    seq.add_annotation(Feature, "CDS", "CDS", [cds])
    seq.add_annotation(Feature, "5'UTR", "5' UTR", [utr])
    return seq


def makeSampleAlignment():
    seq1 = makeSampleSequence()
    seq2 = makeSampleSequence(with_gaps=True)
    seqs = {"FAKE01": seq1, "FAKE02": seq2}
    aln = make_aligned_seqs(data=seqs, array_align=False)
    aln.add_annotation(Feature, "misc_feature", "misc", [(12, 25)])
    aln.add_annotation(Feature, "CDS", "blue", [(15, 25)])
    aln.add_annotation(Feature, "5'UTR", "red", [(2, 4)])
    aln.add_annotation(Feature, "LTR", "fake", [(2, 15)])
    return aln


class TestAnnotations(unittest.TestCase):
    def setUp(self):
        self.seq = makeSampleSequence()
        self.aln = makeSampleAlignment()

    def test_inherit_feature(self):
        """should be able to subclass and extend _Feature"""

        class NewFeat(_Feature):
            def __init__(self, *args, **kwargs):
                super(NewFeat, self).__init__(*args, **kwargs)

            def newMethod(self):
                if len(self.map.spans) > 1:
                    as_one = self.as_one_span()  # should create new instance of NewFeat
                    return as_one.newMethod()
                return True

        seq = DNA.make_seq("ACGTACGTACGT")
        f = seq.add_annotation(
            NewFeat, as_map([(1, 3), (5, 7)], len(seq)), type="gene", name="abcd"
        )
        self.assertEqual(type(f.as_one_span()), NewFeat)
        self.assertEqual(type(f.get_shadow()), NewFeat)
        f2 = seq.add_annotation(
            NewFeat, as_map([(3, 5)], len(seq)), type="gene", name="def"
        )

        self.assertEqual(
            type(seq.get_region_covering_all([f, f2], feature_class=NewFeat)), NewFeat
        )
        # now use the new method
        f.newMethod()

    def test_slice_seq_with_annotations(self):
        newseq = self.seq[:5] + self.seq[10:]
        for annot_type in ["CDS", "5'UTR"]:
            orig = str(list(self.seq.get_by_annotation(annot_type))[0])
            new = str(list(newseq.get_by_annotation(annot_type))[0])
            assert orig == new, (annot_type, orig, new)

    def test_aln_annotations(self):
        """test that annotations to alignment and its' sequences"""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"},
        }
        seq_expecteds = {
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "TTT", "FAKE02": "TTT"},
        }
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = list(self.aln.get_by_annotation(annot_type))[0].to_dict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, (annot_type, expected, observed)
            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it
            for name in self.aln.names:
                observed = list(
                    self.aln.named_seqs[name].data.get_by_annotation(annot_type)
                )[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (annot_type, name, expected, observed)

    def test_slice_aln_with_annotations(self):
        """test that annotations of sequences and alignments survive alignment
        slicing."""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCTTTTT", "FAKE02": "CCCTTTTT"},
        }
        newaln = self.aln[:5] + self.aln[10:]
        feature_list = newaln.get_annotations_matching("LTR")
        for annot_type in ["LTR", "misc_feature", "CDS", "5'UTR"]:
            feature_list = newaln.get_annotations_matching(annot_type)
            new = newaln.get_region_covering_all(feature_list).get_slice().to_dict()
            expected = aln_expecteds[annot_type]
            assert expected == new, (annot_type, expected, new)
            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it
            for name in self.aln.names:
                orig = str(
                    list(self.aln.get_annotations_from_seq(name, annot_type))[
                        0
                    ].get_slice()
                )
                new = str(
                    list(newaln.get_annotations_from_seq(name, annot_type))[
                        0
                    ].get_slice()
                )
                assert orig == new, (name, annot_type, orig, new)

    def test_feature_projection(self):
        expecteds = {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"}
        aln_ltr = self.aln.get_annotations_matching("LTR")[0]
        for seq_name in ["FAKE01", "FAKE02"]:
            expected = expecteds[seq_name]
            seq_ltr = self.aln.project_annotation(seq_name, aln_ltr)
            if "-" in expected:
                self.assertRaises(ValueError, seq_ltr.get_slice)
                seq_ltr = seq_ltr.without_lost_spans()
                expected = expected.replace("-", "")
            self.assertEqual(seq_ltr.get_slice(), expected)

    def test_feature_copy_annotations_to(self):
        """test correct copy of annotations"""
        orig = DnaSequence("TTTTTTTTTTAAAA", name="Orig")
        annot = orig.add_annotation(Feature, "exon", "fred", [(0, 14)])
        seq = RnaSequence("UUUUUUUUUUAAAA", name="Test")
        got = annot.copy_annotations_to(seq)
        self.assertEqual(len(orig.annotations), len(got.annotations))
        for src, dest in zip(orig.annotations, got.annotations):
            self.assertEqual(src.get_coordinates(), dest.get_coordinates())
            self.assertIsInstance(src, dest.__class__)
            self.assertIs(dest.parent, seq)
        with self.assertRaises(AssertionError):
            _ = annot.copy_annotations_to(seq[:-2])

    def test_reverse_complement(self):
        """test correct translation of annotations on reverse complement."""
        aln_expecteds = {
            "misc_feature": {"FAKE01": "TTTGGGGGGGGGG", "FAKE02": "TTTGGGGGGGGGG"},
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "CC", "FAKE02": "CC"},
            "LTR": {"FAKE01": "CCCAAAATTTTTT", "FAKE02": "CCC-----TTTTT"},
        }

        seq_expecteds = {
            "CDS": {"FAKE01": "GGGGGGGGGG", "FAKE02": "GGGGGGGGGG"},
            "5'UTR": {"FAKE01": "TTT", "FAKE02": "TTT"},
        }

        rc = self.aln.rc()
        # rc'ing an Alignment or Sequence rc's their annotations too. This means
        # slicing returns the same sequence as the non-rc'd alignment/seq
        for annot_type in ["misc_feature", "CDS", "5'UTR", "LTR"]:
            observed = list(self.aln.get_by_annotation(annot_type))[0].to_dict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, ("+", annot_type, expected, observed)
            observed = list(rc.get_by_annotation(annot_type))[0].to_dict()
            expected = aln_expecteds[annot_type]
            assert observed == expected, ("-", annot_type, expected, observed)

            if annot_type in ["misc_feature", "LTR"]:
                continue  # because seqs haven't been annotated with it
            for name in self.aln.names:
                observed = list(
                    self.aln.named_seqs[name].data.get_by_annotation(annot_type)
                )[0]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (
                    "+",
                    annot_type,
                    name,
                    expected,
                    observed,
                )
                observed = list(rc.named_seqs[name].data.get_by_annotation(annot_type))[
                    0
                ]
                observed = str(observed)
                expected = seq_expecteds[annot_type][name]
                assert str(observed) == expected, (
                    "-",
                    annot_type,
                    name,
                    expected,
                    observed,
                )


class TestMapSpans(unittest.TestCase):
    """Test attributes of Map & Spans classes critical to annotation
    manipulation."""

    def test_span(self):
        forward = Span(20, 30)
        reverse = Span(70, 80, reverse=True)
        assert forward.reversed_relative_to(100) == reverse
        assert reverse.reversed_relative_to(100) == forward

    def test_map(self):
        """reversing a map with multiple spans should preserve span relative
        order"""
        forward = [Span(20, 30), Span(40, 50)]
        fmap = Map(spans=forward, parent_length=100)
        fmap_reversed = fmap.nucleic_reversed()
        reverse = [Span(70, 80, reverse=True), Span(50, 60, reverse=True)]
        rmap = Map(spans=reverse, parent_length=100)
        for i in range(2):
            self.assertEqual(fmap_reversed.spans[i], rmap.spans[i])


def test_make_gff3_db():
    """Test if the gff3 database is correctly created
    in terms of the names of columns & the number of entries"""

    from cogent3.core.annotation import GffAnnotationDb
    from cogent3.parse.gff import gff_parser

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb(gff_path)

    """Test the name of columns created"""
    sql_columns = []
    data = db.db.execute("""SELECT * FROM GFF""")
    for column in data.description:
        sql_columns.append(column[0])
    got = sql_columns
    parsed_gff = gff_parser(gff_path)
    expected = list((list(parsed_gff)[0]).keys())[:-1]
    assert got == expected


def test_make_sql_query():
    """Test that the SQLlite queries are correctly formed"""

    from cogent3 import load_unaligned_seqs
    from cogent3.core.annotation import GffAnnotationDb
    from cogent3.parse.gff import gff_parser

    fasta_path = "/Users/kiratalreja/Desktop/short.fa"
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb(gff_path)

    """only bio_type provided"""
    expected = (
        "SELECT * FROM GFF WHERE SeqID == ? AND Start >= ? AND End <= ? AND Type == ?",
        [0, 13720, "CDS"],
    )
    got = db._make_sql_query(start=0, end=len(seq), bio_type="CDS")
    assert expected == got

    """only identifier provided"""
    expected = (
        "SELECT * FROM GFF WHERE SeqID == ? AND Start >= ? AND End <= ? AND Attributes like ?",
        [0, 13720, "%RandomID%"],
    )
    got = db._make_sql_query(start=0, end=len(seq), identifier="RandomID")
    assert expected == got

    """both bio_type and identifier provided"""
    expected = (
        "SELECT * FROM GFF WHERE SeqID == ? AND Start >= ? AND End <= ? AND Type == ? AND Attributes like ?",
        [0, 13720, "CDS", "%RandomID%"],
    )
    got = db._make_sql_query(
        start=0, end=len(seq), bio_type="CDS", identifier="RandomID"
    )
    assert expected == got

    """check exception when both bio_type and identifier are missing"""
    import pytest

    with pytest.raises(ValueError):
        db._make_sql_query(start=0, end=len(seq))


def test_populate_from_file():
    """Test that the database is populated with the correct
    number of columns"""
    from cogent3.core.annotation import GffAnnotationDb
    from cogent3.parse.gff import gff_parser

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb()
    db.populate_from_file(gff_path)

    """test the number of rows populated"""
    db.db.execute(""" SELECT * FROM GFF WHERE SeqID=='sequence001' """)
    got = len(list(db.db.fetchall()))
    parsed_gff = gff_parser(gff_path)
    expected = len(list(parsed_gff))
    assert got == expected


def test_db_query():

    """Test that the SQL query returns the correct
    number of rows for different combinations of bio_type/identifier"""

    # PENDING - NEED TO COVER THE CASE WITH NO ROWS RETURNED

    from cogent3 import load_unaligned_seqs
    from cogent3.core.annotation import GffAnnotationDb

    fasta_path = "/Users/kiratalreja/Desktop/short.fa"
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb()
    db.populate_from_file(gff_path)

    """multiple hits for the same identifier"""
    got = len(db.db_query(start=0, end=len(seq), identifier="CDS4"))
    expected = 2
    assert got == expected

    """query for an ID and recieve the children to the ID along with the parent"""
    got = len(db.db_query(start=0, end=len(seq), identifier="gene4"))
    expected = 3
    assert got == expected

    """query for an ID, with no children"""
    got = len(db.db_query(start=0, end=len(seq), identifier="trna1"))
    expected = 1
    assert got == expected

    """query for a bio type, with multiple hits"""
    got = len(db.db_query(start=0, end=len(seq), bio_type="gene"))
    expected = 8
    assert got == expected

    """query for an ID & a bio type, with a single hit"""
    got = len(db.db_query(start=0, end=len(seq), identifier="gene0", bio_type="CDS"))
    expected = 1
    assert got == expected

    """query for an ID & a bio type, with multiple hits"""
    got = len(db.db_query(start=0, end=len(seq), identifier="CDS4", bio_type="CDS"))
    expected = 2
    assert got == expected


def test_find_records():

    """Test that the coordinates the grouped correctly, and features
    formed properly"""

    from cogent3 import load_unaligned_seqs
    from cogent3.core.annotation import GffAnnotationDb

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb()
    db.populate_from_file(gff_path)
    fasta_path = "/Users/kiratalreja/Desktop/short.fa"
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    """combine rows with the same ID"""
    got = len(db.find_records(start=0, end=len(seq), identifier="CDS4"))
    expected = 1
    assert got == expected

    """combine rows with the same ID, when bio_type given"""
    got = len(db.find_records(start=0, end=len(seq), bio_type="CDS"))
    expected = 6
    assert got == expected

    """combine rows with the same ID, when children rows are fetched with the parent"""
    got = len(db.find_records(start=0, end=len(seq), identifier="gene4"))
    expected = 2
    assert got == expected

    """unique ID, single row returned"""
    got = len(db.find_records(start=0, end=len(seq), identifier="id020000"))
    expected = 1
    assert got == expected


if __name__ == "__main__":
    unittest.main()
