# Architecture of the repository

In `data`you will find all the sample files for trexio computations.
In `lib`you will find the executable.
In `src`you will find the source files for the project.

# Sources

The main source is `main.c`. It contains the main function and private functions to compute different things using trexio (repulsion energy, orbitals, ...). A documentation for those functions is available at the beginning of the file.

# How to run

To compile and run, you just need to launch the `compile_and_run.sh` file. To use it, you need to call it with two arguments:
`compile_and_run.sh <file_path> <molecule>` where `file_path` is the TREXIO file to analyse and `molecule` is the name of the molecule analysed, used for printing purposes.