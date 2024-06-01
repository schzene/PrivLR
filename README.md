# PrivLR: Privacy Preserving Logistic Regression

PrivLR is Privacy Preserving Logistic Regression implement with C++, using paillier.

# Build

```bash
mkdir build && cd build
cmake ../
make
```

if you want to turn off the test sample, run:

```bash
mkdir build && cd build
cmake ../ -DPrivLR_TEST=OFF
make
```

# Acknowledgments

This repository includes code from the following external repositories:

[emp-toolkit/emp-tool](https://github.com/emp-toolkit/emp-tool) for Network IO