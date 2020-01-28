Release 0.1-0 (1/28/2020):
New:
  * Created spmat, spvec, and dvec classes
  * Added to spar::reduce namespace:
    - dense() for sparse matrix (all)reduce with dense column intermediary.
    - gather() for sparse matrix (all)reduce with an MPI gather strategy.
  * Added to spar::conv namespace:
    - Eigen converters:
      - spmat_to_eigen()
      - eigen_to_spmat()
      - eigen_to_dvec()
    - R dgCMatrix converters:
      - s4col_to_spvec()
      - spmat_to_s4()
  * Added to spar::gen namespace:
    - bandish() for random sparse matrices that are "nearly" banded
API Changes: None
Bug Fixes: None
Documentation:
  * Added documentation for all non-internal classes/methods except.
