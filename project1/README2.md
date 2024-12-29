# Introduction

TO BE COMPLETEDThis program reads molecular data from a '.h5' file and
computes the Hartree-Fock (HF) energy and Møller--Plesset second-order
perturbation (MP2) energy.

# Theory

The HF energy is given by:
<<<<<<< HEAD
\[E_{\text{HF}} = E_{NN} + 2\sum_{i=1}^{N_{\text{occ}}}\langle i|\hat{h}|i\rangle + \sum_{i=1}^{N_{\text{occ}}}\sum_{j=1}^{N_{\text{occ}}} \left(2\langle ij|ij\rangle - \langle ij|ji\rangle\right)\]

The MP2 energy is then calculated as:
\[E_{\text{MP2}} = E_{\text{HF}} + 
\sum_{(i,j)\in\text{occupied}}\sum_{(a,b)\in\text{virtual}}\langle ij|ab\rangle\frac{2\langle ij|ab\rangle-\langle ij|ba\rangle}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}\]
=======

![E\_{\\text{HF}} = E\_{NN} + 2\\sum\_{i=1}\^{N\_{\\text{occ}}}\\langle i\|\\hat{h}\|i\\rangle + \\sum\_{i=1}\^{N\_{\\text{occ}}}\\sum\_{j=1}\^{N\_{\\text{occ}}} \\left(2\\langle ij\|ij\\rangle - \\langle ij\|ji\\rangle\\right)](https://latex.codecogs.com/png.latex?E_%7B%5Ctext%7BHF%7D%7D%20%3D%20E_%7BNN%7D%20%2B%202%5Csum_%7Bi%3D1%7D%5E%7BN_%7B%5Ctext%7Bocc%7D%7D%7D%5Clangle%20i%7C%5Chat%7Bh%7D%7Ci%5Crangle%20%2B%20%5Csum_%7Bi%3D1%7D%5E%7BN_%7B%5Ctext%7Bocc%7D%7D%7D%5Csum_%7Bj%3D1%7D%5E%7BN_%7B%5Ctext%7Bocc%7D%7D%7D%20%5Cleft%282%5Clangle%20ij%7Cij%5Crangle%20-%20%5Clangle%20ij%7Cji%5Crangle%5Cright%29 "E_{\text{HF}} = E_{NN} + 2\sum_{i=1}^{N_{\text{occ}}}\langle i|\hat{h}|i\rangle + \sum_{i=1}^{N_{\text{occ}}}\sum_{j=1}^{N_{\text{occ}}} \left(2\langle ij|ij\rangle - \langle ij|ji\rangle\right)")

The MP2 energy is then calculated as:

![E\_{\\text{MP2}} = E\_{\\text{HF}} + 
\\sum\_{(i,j)\\in\\text{occupied}}\\sum\_{(a,b)\\in\\text{virtual}}\\langle ij\|ab\\rangle\\frac{2\\langle ij\|ab\\rangle-\\langle ij\|ba\\rangle}{\\varepsilon\_i+\\varepsilon\_j-\\varepsilon\_a-\\varepsilon\_b}](https://latex.codecogs.com/png.latex?E_%7B%5Ctext%7BMP2%7D%7D%20%3D%20E_%7B%5Ctext%7BHF%7D%7D%20%2B%20%0A%5Csum_%7B%28i%2Cj%29%5Cin%5Ctext%7Boccupied%7D%7D%5Csum_%7B%28a%2Cb%29%5Cin%5Ctext%7Bvirtual%7D%7D%5Clangle%20ij%7Cab%5Crangle%5Cfrac%7B2%5Clangle%20ij%7Cab%5Crangle-%5Clangle%20ij%7Cba%5Crangle%7D%7B%5Cvarepsilon_i%2B%5Cvarepsilon_j-%5Cvarepsilon_a-%5Cvarepsilon_b%7D "E_{\text{MP2}} = E_{\text{HF}} + 
\sum_{(i,j)\in\text{occupied}}\sum_{(a,b)\in\text{virtual}}\langle ij|ab\rangle\frac{2\langle ij|ab\rangle-\langle ij|ba\rangle}{\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b}")
>>>>>>> lav

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
