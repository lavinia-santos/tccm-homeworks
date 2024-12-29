# Introduction

TO BE COMPLETEDThis program reads molecular data from a ‘.h5‘ file and
computes the Hartree-Fock (HF) energy and Møller–Plesset second-order
perturbation (MP2) energy.

# Theory

The HF energy is given by:
\[E_{\text{HF}} = E_{NN} + 2\sum_{i=1}^{N_{\text{occ}}}\langle i|\hat{h}|i\rangle + \sum_{i=1}^{N_{\text{occ}}}\sum_{j=1}^{N_{\text{occ}}} \left(2\langle ij|ij\rangle - \langle ij|ji\rangle\right)\]

The MP2 energy is then calculated as:
\[E_{\text{MP2}} = E_{\text{HF}} + 
\sum_{(i,j)\in\text{occupied}}\sum_{(a,b)\in\text{virtual}}\langle ij|ab\rangle\frac{2\langle ij|ab\rangle-\langle ij|ba\rangle}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}\]

# Dependencies

The program depends on:

  - **HDF5 library**: Install using `sudo apt-get install libhdf5-dev`.

  - **TREXIO library**: Install following the instructions on [TREXIO
    GitHub](https://github.com/TREX-CoE/trexio).

# Installation and Usage

Clone the repository, navigate to the project directory, and compile
using the following commands:

    $ git clone git@github.com:yourusername/yourrepo.git
    $ cd yourrepo
    $ ./compile.sh

Run the program as follows:

    ./main.o [path/to/input.h5] [molecule_name]

# Example

    $ ./main.o ./data/h2o.h5 H2O

Expected output:

    HF energy: -76.026799
    MP2 energy: -76.230759

# Testing

Use the provided `test.sh` script to validate the program:

    $ ./test.sh
