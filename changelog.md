# Changes since release 2022.8.24a1

## Contributors

Thanks to our contributors!

### Accepted PRs from

- Gavin Huttley
- KatherineCaley
- Nick Shahmaras
- Xingjian Leng

### Identified a Bug

- StephenRogers1

## ENH

Significant refactor of composable apps. This is the backwards incompatible change we warned of in the last release! We now use a decorator `define_app` (on classes or functions) instead of class inheritance. Please see [the c3dev wiki](https://github.com/cogent3/cogent3/wiki/composable-functions) for examples on how to port from old-style to new-style composable apps.

We updated to the latest NCBI versions of genetic codes. Note, the name of genetic code 1 has changed from "Standard Nuclear" to "Standard".

## BUG

- Fix progressive alignment bug when a guide-tree with zero edge lengths was encountered.
- Non-stationary independent tuple models can now be serialised.

## DEP

- We have removed support for python 3.7.
- We have made scipy a dependency and begun deprecating statistical functions that are available in scipy. All deprecated functions have a warning that indicates the scipy replacement. The deprecated functions are: combinations, chi_high, chdtri, z_high, z_low function, chi_low, binomial_high, binomial_low, f_high, f_low, t_low and t_high.

# Changes since release 2022.5.25a1

## Contributors

- Gavin Huttley
- Nick Shahmaras

## Notice of upcoming major changes

The definition of composable apps by inheritance from the `Composable` app base class will no longer be supported from 2022.11. This is to be replaced by a decorator which greatly simplifies constructing new composable apps. Guidance on how to port existing code will be [posted to a dedicated cogent3 wiki page](https://github.com/cogent3/cogent3/wiki/composable-functions). The `cogent3` composable apps will continue to be available and to work as per usual.

We will drop support for python 3.7 by release 2022.10.

## API

- The following method names are marked for deprecation. `Sequence.gettype` for `Sequence.get_type`, `Sequence.resolveambiguities` for `Sequence.resolved_ambiguities`.

## ENH

- Refactor dotplot, fixes #1060. This is a major performance improvement in dotplot code, making it much faster for long DNA sequences. I've implemented a seed-and-extend algorithm which identifies matching k-mers between sequences and only works on extending those. The key interfaces are `find_matched_paths()` which takes a `SeqKmers()` instance, the two sequences to be compared and returns a `MatchedSeqPaths()` instance.
- Added Sequence methods to compute all k-mers, fixes #1012
- `Sequence.iter_kmers()`, generator yielding all overlapping k-mers
- Enhanced placement of dendrogram scale bars
- `cogent3.open_` now handles urls too! This enables files to be downloaded from the web using the convenience functions, e.g. `load_aligned_seqs("https://path/to/somefile.fasta")` will now work. Many thanks to Nick Shahmaras for assistance in doing this!

# Changes since release 2022.4.20a1

## Contributors

- Gavin Huttley

## ENH

- new `cogent3.util.parallel.as_completed()` generator function
    - `as_completed()` wraps MPI or `concurrent.futures` executors and delivers results as they are completed. In contrast, `parallel.imap()` / `parallel.map()` deliver results in the same order as the input series. The advantage of `as_completed()` is the interval of result arrival at the parent process is better distributed.
- new function `cogent3.load_seq()` loads a single sequence from a file
- convert substitution model `__str__` to `__repr__`; more useful since `__repr__` is called also by str().

## BUG

- fixes to `annotation_from_gff()` method on annotatable sequence / alignment objects
    - method would break if GFF records had no ID. This situation is quite common in some Ensembl gff3 files. We generate a "no-id-#" identifier in those cases.
    - we now add features are added to their parent feature.
- improve consistency in setting motif_probs on likelihood function
    - only apply a pseudocount if optimising motif probs and at least one state has zero frequency, default pseudocount is 0.5. Thanks to StephenRogers1 for finding this issue!

## DOC

- document the API of the new `load_seq()` function

# Changes since release 2022.4.15a1

## DEP

- added warning that we will drop support for python 3.7 by 2022.10. This means the developer version will switch to python 3.8 from 2022.6.
- discontinued delimiter argument from parse.table.load_delimited

# Changes since release 2021.10.12a1

## Contributors

- Gavin Huttley
- u6675275

## API

- moved all io related functions classes from util.misc to util.io, indicating their removal after version 2022.4
- app.result objects require source instance of str or pathlib.Path
- fail if users set motif prob optimisation via sm_args in app.evo.model as value is over ridden by the explicit argument, need to block this as effect is major

## BUG

- RichGenbankParser moltype argument now overrides file spec, if provided, this defines the moltype of the returned sequence, otherwise the moltype is determined from the genbank file meta-data
- fix initialise from nested params for codon models
- load_tree now handles pathlib.Path's as input, fixes #991
- writer composable apps apply_to now handles provided logger
- fixed serialisation of multi-locus likelihood functions with constrained motif probs
- support multiple calls of to_rich_dict()
- solve case where optimiser gets an invalid starting vector
- solved case where optimised parameter values are outside bounds

## DEP

- removed deprecated function for median, use numpy.median instead
- removed deprecated index argument from table constructors, use index_name instead
- cogent3.math periodicity classes method names pep8, old names retained, with deprecation warnings

## ENH

- Drawable.plotly_figure property returns plotly graph object Figure instance
- refactor of cogent3.app.composable.appify so decorated functions can be pickled
- app.evo.model handles sequential fitting of models with mixed process. Sequential fitting now works if lf_args includes specifying edges for using a discrete-time Markov process
- add optimise_motif_probs argument to app.evo.model
- add upper argument to app.evo.model
- now support python 3.10
- added register_model decorator to cogent3.evolve.models. Used for simplifying discovery of canned substitution models. Users can now use this mechanism too for adding their own custom models. Doing this smoothes usage of custom models with cogent3.app.evo.model. A further benefit is the inclusion of a model to the appropriate module attributes is now done automatically.
- generalise Jensen-Shannon calculations to > 2 distributions
- the register_deserialiser class takes a series of strings that serve to uniquely identify the "type" value in a dict to be reconstituted using the decorated function. This enables support for user defined custom json storage.
- add type hint for input paths to most commonly used loaders
- time-heterogeneity support mixed discrete and continuous-time models
- more compact representation of datastore summary_incomplete
- more refinements on summary_logs
- cogent3.app.io.register_datastore_reader enables development of third party readers / loaders to be developed. Registering a reader class requires decorating it with the filename suffix that will distinguish that content type. Still limited to reading from files only.
- improve general stationary model numerical precision tolerance

# Since release 2021.5.7a1

## Contributors

- GavinHuttley

## DEV

- added missing `dev` requires-extras to pyproject.toml for installing all packages required for development

# Since release 2021.5.7a

## Contributors

- GavinHuttley
- jamesmartini
- KatherineCaley

## API

- ValueError if any tips missing in TreeNode.lowest_common_ancestor()
- added index_name argument to Table.to_categorical(), allows specifying the category column and getting the categorical table in one statement.

## BUG

- DataStore.write() requires identifiers end with indicated suffix
- cogent3.app.tree.quicktree() now works for 2 sequences
- Alignment.degap() now preserves sequence names
- cogent3.app.io.load_aligned() handles paml format
- fast_slow_dist results can now be saved by write_tabular, a DistanceMatrix.source attribute is created on-the-fly by the fast_slow_dist calculator, enabling it be written
- Alignment.variable_positions(), always report a position as variable if > 1 non-gap characters are present
- SequenceCollection.dotplot() method defaults handle single sequence
- app.evo.boostrap() can now be composed,

## DEV

- change to using flit for package management. This change requires you `python -m pip install flit`. If you clone this repository and want to do a developer install, you should first remove your existing one

    ````bash
    $ python -m pip uninstall cogent3
    ````

    then

    ```bash
    $ flit install -s --python `which python`
    ```

## DEP

- removed WritableZippedDataStore, the zip archive format is inefficient for incremental inclusion of files. Use a tinydb instead.
- replaced interleave_len argument with wrap in sequence format writers
- removed Table.to_rich_html() method, use Table.to_html() instead
 
## ENH

- More robust alignment to reference algorithm. Builds a multiple sequence alignment from a series of pairwise alignments to a reference sequence. cogent3.app.align.align_to_ref() now retains gaps in the reference. This will be modestly slower than previously, but avoids losing information if the choice of reference sequence is a bad one.
- cogent3.app.composable.appify decorator class, simplifies converting a user defined function into a cogent3 composable app
- JSD calculation now uses more accurate math.fsum()