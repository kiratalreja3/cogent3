from memory_profiler import profile

from cogent3 import load_unaligned_seqs
from cogent3.core.annotation import GffAnnotationDb


@profile
def new_approach():
    gff_path = "/Users/kiratalreja/Downloads/Homo_sapiens.GRCh38.108.chromosome.1.gff3"
    db = GffAnnotationDb(gff_path)

    fasta_path = "/Users/kiratalreja/Downloads/Homo_sapiens.GRCh38.dna.chromosome.1.fa"
    seqs = load_unaligned_seqs(fasta_path, moltype="dna")
    seq = seqs.seqs[0]
    seq.name = "1"

    seq.annotate_from_db(db, bio_type="gene")
    seq.get_annotations_matching("gene")


if __name__ == "__main__":
    new_approach()
