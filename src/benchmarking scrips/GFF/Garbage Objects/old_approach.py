import gc

from cogent3 import load_unaligned_seqs
from cogent3.parse.gff import gff_parser


def old_approach():
    def fix_gff(gff):
        return [r for r in gff if r["Attributes"]["ID"]]

    fasta_path = "/Users/kiratalreja/Downloads/Homo_sapiens.GRCh38.dna.chromosome.1.fa"
    gff3_path = "/Users/kiratalreja/Downloads/Homo_sapiens.GRCh38.108.chromosome.1.gff3"

    gff = list(gff_parser(gff3_path))
    gff = fix_gff(gff)

    seqs = load_unaligned_seqs(fasta_path, moltype="dna")

    seq = seqs.seqs[0]
    seq.name = "1"

    seq.annotate_from_gff(gff, pre_parsed=True)
    seq.get_annotations_matching("gene")


if __name__ == "__main__":
    old_approach()
    print(len(gc.get_objects()))
