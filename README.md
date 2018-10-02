## Estimating the density of states using nonequilibrium trajectories.

Grant M. Rotskoff & Eric Vanden-Eijnden: [arXiv:1809.11132](https://arxiv.org/abs/1809.11132)

This code can be used to perform density of states calculations, as described in the arXiv reference above.

### Build instructions

Dependencies: armadillo, hdf5, C++-11

Compiling is simple, e.g.,

`make quartic`

Usage is specified within the run scripts, but is typically

`./run [parameters]`


Additional potentials can be implemented with header files like `quartic.h`
