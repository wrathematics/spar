# Release 0.2-0 (//):

New:
  * Added to spar::gen namespace
    - rand() for uniform random sparse matrices
    - banded() for banded sparse matrices
  * Added two benchmarks to the `benchmarks/` tree.

API Changes:
  * All headers are now contained in the `spar/` tree.
      - `src/spar.hpp` becomes `src/spar/spar.hpp`
      - `src/reduce.hpp` becomes `src/spar/reduce.hpp`
      - etc.
  * All classes are now in the spar namespace.
      - `spmat<int, float>` becomes `spar::spmat<int, float>`
      - `spvec<int, float>` becomes `spar::spvec<int, float>`
      - etc.

Bug Fixes:
  * Fixed sparsity/density lookup (was backwards).
  * Fixed some edge case memory errors with the initial sparse matrix size in
    the reducers.
  * Initial sparse matrix memory size is less aggressive in reducers.

Documentation:
  * Fixed documentation generation for spar::gen.
  * Fixed some minor discrepancies in spvec documentation.





# Release 0.1-0 (1/28/2020):

New:
  * Created spmat, spvec, and dvec classes
  * Created spar::reduce namespace
  * Added to spar::reduce namespace
    - dense() for sparse matrix (all)reduce with dense column intermediary.
    - gather() for sparse matrix (all)reduce with an MPI gather strategy.
  * Created spar::conv namespace
  * Added to spar::conv namespace
    - Eigen converters:
      - spmat_to_eigen()
      - eigen_to_spmat()
      - eigen_to_dvec()
    - R dgCMatrix converters:
      - s4col_to_spvec()
      - spmat_to_s4()
  * Created spar::gen namespace
  * Added to spar::gen namespace
    - bandish() for random sparse matrices that are "nearly" banded

API Changes: None

Bug Fixes: None

Documentation:
  * Added documentation for all non-internal classes/methods except.
