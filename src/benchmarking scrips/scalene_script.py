from pycallgraph2 import PyCallGraph
from pycallgraph2 import Config
from pycallgraph2 import GlobbingFilter
from pycallgraph2.output import GraphvizOutput
graphviz = GraphvizOutput()
graphviz.output_file = 'basic1.png'
with PyCallGraph(output=graphviz):
        from cogent3.core.annotation import GffAnnotationDb
        from cogent3.parse.gff import gff_parser

        gff_path = '/Users/kiratalreja/Downloads/Homo_sapiens.GRCh38.108.chromosome.1.gff3'
        db = GffAnnotationDb(gff_path)
        db.find_records(bio_type='gene')