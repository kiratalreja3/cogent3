import unittest


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2007-2022, The Cogent Project"
__credits__ = ["Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "2022.8.24a1"
__maintainer__ = "Gavin Huttley"
__email__ = "gavin.huttley@anu.edu.au"
__status__ = "Production"


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


def test_ordered_values():
    from cogent3.core.annotation_db import _ordered_values
    parsed_gff = gff_parser(gff_path)
    example_dict = list(parsed_gff)[0]  
    got = _ordered_values(example_dict)
    expected = ['sequence001', 'mine','gene', 189, 255, '.', '+', '.', {'ID': 'gene0','Dbxref': 'ASAP:ABE-0000006','gene': 'thrL','gene_synonym': 'ECK0001'}]
    assert got == expected

    def return_type(values):
        types = []
        for v in values:
            types.append(type(v))
        return types

    got = return_type(expected)
    expected = [str, str, str, int, int, str, str, str, dict]
    assert got == expected


if __name__ == "__main__":
    test_db_query()
    test_find_records()
    test_make_gff3_db()
    test_make_sql_query()
    test_populate_from_file()
    test_ordered_values()