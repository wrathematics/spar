# spvec

* **Version:** 0.1-0
* **Status:** [![Build Status](https://travis-ci.org/wrathematics/spvec.png)](https://travis-ci.org/wrathematics/spvec)
* **License:** [BSL-1.0](http://opensource.org/licenses/BSL-1.0)
* **Project home**: https://github.com/wrathematics/spvec
* **Bug reports**: https://github.com/wrathematics/spvec/issues


spvec is a small, header-only C++ library with a very specific purpose. You very likely will not find this useful; I'm only putting it here in its own repository because I need lots of tests that I don't want cluttering my main project.

The purpose of the library to allow one to add many sparse matrices in CSC format, with the matrices spread across multiple processors. If you do not need to do this very specific thing, then you should use something else like [Eigen](http://eigen.tuxfamily.org/) or [Armadillo](http://arma.sourceforge.net/) for your sparse matrices.



## Dependencies and Tests

There are no external dependencies. Tests use [catch2](https://github.com/catchorg/Catch2), a copy of which is included under `tests/catch`.

To build the tests, modify `tests/make.inc` as appropriate and type `make`.
