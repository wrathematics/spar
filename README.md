# spar

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/spar.png)](https://travis-ci.org/wrathematics/spar)
* **License:** [BSL-1.0](http://opensource.org/licenses/BSL-1.0)
* **Project home**: https://github.com/wrathematics/spar
* **Bug reports**: https://github.com/wrathematics/spar/issues
* **Documentation**: http://librestats.com/spar/html/index.html

<img align="right" src="./docs/logo/spar_med.png" />

spar is the Sparse Allreduce library, a header-only C++14 framework. Its purpose is to enable the addition of many sparse matrices in CSC format, with the matrices spread across multiple processors.

A sparse allreduce is different from the addition of two distributed, sparse matrices. In that case, every process contains a piece of each of the two distributed matrices. After the sum each process contains a piece of the new distributed matrix, whose entries are the sum of the two. In a sparse allreduce, each of the, say, `p` processes has a sparse matrix, all of the `p` matrices are added, and afterwards, one or all processes contain the sum.



## Dependencies and Tests

If you want to use the actual reducers, you will need an installation of MPI.

Tests use [catch2](https://github.com/catchorg/Catch2), a copy of which is included under `tests/catch`. To build the tests, modify `tests/make.inc` as appropriate and type `make`.

To use the library with Eigen, you will need a copy of [the Eigen headers](http://eigen.tuxfamily.org/index.php?title=Main_Page). You may find the [eigen example](examples/mpi/eigen.cpp) helpful.
