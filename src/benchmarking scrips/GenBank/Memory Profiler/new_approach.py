from memory_profiler import profile

@profile
def new_approach():
    from cogent3.core.annotation_db import GenbankAnnotationDb
    gff_path = "/Users/kiratalreja/Downloads/NC_000913.3.gb"
    db = GenbankAnnotationDb(gff_path)
    db.find_records(bio_type="gene")


if __name__ == "__main__":
    new_approach()
