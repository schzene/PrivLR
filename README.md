# PrivLR: Privacy Preserving Logistic Regression

PrivLR is Privacy Preserving Logistic Regression implement with C++, using paillier.

# Build

The following library must be installed before you build PrivLR:

* GMP
* NTL

You can build PrivLR by following command: 

```bash
mkdir build && cd build
cmake ../
make
```

If you want to turn off the test sample, run:

```bash
mkdir build && cd build
cmake ../ -DPrivLR_TEST=OFF
make
```

# Acknowledgments

This repository includes code from the following external repositories:

[mpc-msri/EzPC](https://github.com/mpc-msri/EzPC) for Network IO