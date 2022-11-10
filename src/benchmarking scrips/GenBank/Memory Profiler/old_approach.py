

@profile
def old_approach():
    from cogent3 import load_unaligned_seqs
    genbank_path = "/Users/kiratalreja/Downloads/NC_000913.3.gb"
    seqs = load_unaligned_seqs(genbank_path, moltype="dna")
    
if __name__ == "__main__":
    old_approach()
