from cogent3.core.annotation import GffAnnotationDb


gff_path = '/Users/kiratalreja/Downloads/Homo_sapiens.GRCh38.108.chromosome.1.gff3'
db = GffAnnotationDb(gff_path)
db.find_records(bio_type='gene')

