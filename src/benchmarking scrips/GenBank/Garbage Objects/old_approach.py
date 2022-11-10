import gc

from cogent3 import load_unaligned_seqs
from cogent3.parse.gff import gff_parser


def old_approach():
    from cogent3 import load_unaligned_seqs
    genbank_path = "/Users/kiratalreja/Downloads/NC_000913.3.gb"
    seqs = load_unaligned_seqs(genbank_path, moltype="dna")
    seq = seqs.seqs[0]
    #seq.get_annotations_matching("gene")

if __name__ == "__main__":
    old_approach()
    print(len(gc.get_objects()))
