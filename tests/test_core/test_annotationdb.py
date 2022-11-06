from cogent3.core.annotation_db import GffAnnotationDb
from cogent3.core.annotation_db import GenbankAnnotationDb
from cogent3.parse.gff import gff_parser
from cogent3 import open_
from cogent3.parse.genbank import MinimalGenbankParser
from cogent3 import load_unaligned_seqs
gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
genbank_path = "/Users/kiratalreja/Downloads/NC_000913.3.gb"

def test_make_db():
    """Test if the gff3 database is correctly created
    in terms of the names of columns & the number of entries"""

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
    
    sql_columns = []
    db = GenbankAnnotationDb(genbank_path)
    data = db.db.execute("""SELECT * FROM GENBANK""")
    for column in data.description:
        sql_columns.append(column[0])
    got = sql_columns
    expected = ["LocusID","Type","Spans","Locus_Tag","Start","End","Strand"]
    assert got == expected
    
   
def test_make_sql_query():
    """Test that the SQLlite queries are correctly formed"""
    
    fasta_path = "/Users/kiratalreja/Desktop/short.fa"
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb(gff_path)

    """only the bio_type is provided"""
    expected = ('SELECT * FROM GFF WHERE Type == ?', ['gene'])
    got = db._make_sql_query(seq_name=None,bio_type="gene",identifier=None,start=None,end=None)
    assert got == expected 

    """only the identifier is provided"""
    expected = ('SELECT * FROM GFF WHERE Attributes like ?', ['%RandomIdentifier%'])
    got = db._make_sql_query(seq_name=None,bio_type=None,identifier="RandomIdentifier",start=None,end=None)
    assert got == expected

    """both identifier and bio_type provided"""
    expected = ('SELECT * FROM GFF WHERE Type == ? AND Attributes like ?',['CDS', '%RandomIdentifier%'])
    got = db._make_sql_query(seq_name=None,bio_type="CDS",identifier="RandomIdentifier",start=None,end=None)
    assert got == expected

    """start provided along with identifier and bio_type"""
    expected = ('SELECT * FROM GFF WHERE Type == ? AND Attributes like ? AND Start >= ?',['CDS', '%RandomIdentifier%', 0])
    got = db._make_sql_query(seq_name=None,bio_type="CDS",identifier="RandomIdentifier",start=0,end=None)
    assert got == expected

    """end provided along with identifier and bio_type"""
    expected = ('SELECT * FROM GFF WHERE Type == ? AND Attributes like ? AND End < ?',['CDS', '%RandomIdentifier%', 5000])
    got = db._make_sql_query(seq_name=None,bio_type="CDS",identifier="RandomIdentifier",start=None,end=5000)
    assert got == expected

    """start and end provided along with identifier and bio_type"""
    expected = ('SELECT * FROM GFF WHERE Type == ? AND Attributes like ? AND Start >= ? AND End < ?',['CDS', '%RandomIdentifier%', 0, 5000])
    got = db._make_sql_query(seq_name=None,bio_type="CDS",identifier="RandomIdentifier",start=0,end=5000)
    assert got == expected

    """all five attributes provided"""
    expected = ('SELECT * FROM GFF WHERE SeqID == ? AND Type == ? AND Attributes like ? AND Start >= ? AND End < ?', ['1', 'CDS', '%RandomIdentifier%', 0, 5000])
    got = db._make_sql_query(seq_name="1",bio_type="CDS",identifier="RandomIdentifier",start=0,end=5000)
    assert got == expected

    """check exception when both bio_type and identifier are missing"""
    import pytest

    with pytest.raises(ValueError):
        db._make_sql_query(seq_name=None,bio_type=None,identifier=None,start=None,end=None)

    """check exception when both bio_type and identifier are missing, even if other attributes"""

    with pytest.raises(ValueError):
        db._make_sql_query(seq_name="1",bio_type=None,identifier=None,start=0,end=1000)

    """check exception when only seq_name is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(seq_name="1",bio_type=None,identifier=None,start=None,end=None)

    """check exception when only start is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(seq_name=None,bio_type=None,identifier=None,start=0,end=None)
    
    """check exception when only end is provided"""

    with pytest.raises(ValueError):
        db._make_sql_query(seq_name=None,bio_type=None,identifier=None,start=None,end=1000)

def test_populate_from_file():
    """Test that the database is populated with the correct
    number of columns"""
    from cogent3.parse.gff import gff_parser

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb(gff_path)

    """test the number of rows populated"""
    db.db.execute(""" SELECT * FROM GFF """)
    got = len(list(db.db.fetchall()))
    parsed_gff = gff_parser(gff_path)
    expected = len(list(parsed_gff))
    assert got == expected

def test_db_query():

    """Test that the SQL query returns the correct
    number of rows for different combinations of bio_type/identifier"""

    from cogent3 import load_unaligned_seqs

    fasta_path = "/Users/kiratalreja/Desktop/short.fa"
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb(gff_path)

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

    gff_path = "/Users/kiratalreja/Downloads/prok_NoLocusTags.gff"
    db = GffAnnotationDb(gff_path)
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

def test_distinct():

    """Test that the number of unique values returned are correct"""

    db = GffAnnotationDb(gff_path)

    """An empty set when nothing is provided"""
    got = db.distinct()
    expected = set()
    assert got == expected

    """Only the bio_type is provided"""
    got = len(db.distinct(bio_type=True)['Type'])
    expected = 9
    assert got == expected 
    
    """Only the seq_name is provided"""
    got = len(db.distinct(seq_name=True)['SeqID'])
    expected = 1
    assert got == expected
    
    """Only the identifier is provided"""
    got = len(db.distinct(identifier=True)['identifier'])
    expected = 21
    assert got == expected 
    
    """Check that the set has three keys when all three arguments are True"""
    got = len(db.distinct(bio_type=True,seq_name=True,identifier=True))
    expected = 3
    assert got == expected
    
    """Correct number of Type values when all three arguments are True"""
    got = len(db.distinct(bio_type=True,seq_name=True,identifier=True)['Type'])
    expected = 9
    assert got == expected
    
    """Correct number of identifier values when all three arguments are True"""
    got = len(db.distinct(bio_type=True,seq_name=True,identifier=True)['identifier'])
    expected = 21
    assert got == expected
    
    """Correct number of SeqID values when all three arguments are True"""
    got = len(db.distinct(bio_type=True,seq_name=True,identifier=True)['SeqID'])
    expected = 1
    assert got == expected

def test_ordered_values():
    """Check if the dictonary of GFF records is correctly formed into a list for loading into the in-memory database"""
    from cogent3.core.annotation_db import _ordered_values
    
    """Check the values explicitly"""
    parsed_gff = gff_parser(gff_path)
    example_dict = list(parsed_gff)[0]  
    got = _ordered_values(example_dict)
    expected = ['sequence001', 'mine','gene', 189, 255, '.', '+', '.', {'ID': 'gene0','Dbxref': 'ASAP:ABE-0000006','gene': 'thrL','gene_synonym': 'ECK0001'}]
    assert got == expected
    
    """Check the data types of the resulting list"""
    def return_type(values):
        types = []
        for v in values:
            types.append(type(v))
        return types

    got = return_type(_ordered_values(example_dict))
    expected = [str, str, str, int, int, str, str, str, dict]
    assert got == expected



if __name__ == "__main__":
    test_db_query()
    test_find_records()
    test_make_sql_query()
    test_populate_from_file()
    test_distinct()

