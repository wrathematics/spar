# spar

* **Version:** 0.2-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/spar.png)](https://travis-ci.org/wrathematics/spar)
* **License:** [BSL-1.0](http://opensource.org/licenses/BSL-1.0)
* **Project home**: https://github.com/wrathematics/spar
* **Bug reports**: https://github.com/wrathematics/spar/issues
* **Documentation**: http://librestats.com/spar/html/index.html

<img align="right" src="./docs/logo/spar_med.png" />

spar is the Sparse Allreduce library, a header-only C++14 framework. Its purpose is to enable the addition of many sparse matrices in CSC format, with the matrices spread across multiple processors.

A sparse allreduce is different from the addition of two distributed, sparse matrices. In that case, every process contains a piece of each of the two distributed matrices. After the sum each process contains a piece of the new distributed matrix, whose entries are the sum of the two. In a sparse allreduce, each of the, say, `p` processes has a sparse matrix, all of the `p` matrices are added, and afterwards, one or all processes contain the sum.



## Installation and Other Software

The library is header-only so no installation is strictly necessary. You can just include a copy/submodule in your project. However, if you want some analogue of `make install`, then you could do something like:

```bash
ln -s ./src/spar /usr/include/
```

Stable releases are posted [on GitHub](https://github.com/wrathematics/spar/releases). Additionally, you can download the development version via:

```
git clone --recurse-submodules https://github.com/wrathematics/spar.git
```

To use the reducers, you will need an installation of MPI.

Tests use [catch2](https://github.com/catchorg/Catch2), a copy of which is included under `tests/catch`. To build the tests, modify `tests/make.inc` as appropriate and type `make`.

To use the library with Eigen, you will need a copy of [the Eigen headers](http://eigen.tuxfamily.org/index.php?title=Main_Page). You may find the [eigen example](examples/mpi/eigen.cpp) helpful.
