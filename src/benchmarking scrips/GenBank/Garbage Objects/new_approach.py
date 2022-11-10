import gc

from cogent3 import load_unaligned_seqs
from cogent3.core.annotation import GenbankAnnotationDb


def new_approach():
    gff_path = "/Users/kiratalreja/Downloads/NC_000913.3.gb"
    db = GenbankAnnotationDb(gff_path)
    
   


if __name__ == "__main__":
    new_approach()
    print(len(gc.get_objects()))
