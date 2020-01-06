# spar

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/spar.png)](https://travis-ci.org/wrathematics/spar)
* **License:** [BSL-1.0](http://opensource.org/licenses/BSL-1.0)
* **Project home**: https://github.com/wrathematics/spar
* **Bug reports**: https://github.com/wrathematics/spar/issues


spar is a small, header-only C++14 library. The purpose of the library to allow one to add many sparse matrices in CSC format, with the matrices spread across multiple processors.



## Dependencies and Tests

There are no external dependencies. Tests use [catch2](https://github.com/catchorg/Catch2), a copy of which is included under `tests/catch`.

To build the tests, modify `tests/make.inc` as appropriate and type `make`.
