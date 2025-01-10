# Project Overview

This project is designed to read `.h5` files using the TREXIO library to compute the MP2 energy.

## Repository Architecture

- **`data/`**: Contains all the sample files for TREXIO computations.
- **`lib/`**: Contains the compiled executable.
- **`src/`**: Contains the source files for the project, including `main.c`.

## Source Files

The main source file is **`main.c`**, which includes:
- The main function.
- Private functions to compute different properties using TREXIO, such as:
  - Repulsion energy.
  - Molecular orbitals.
  - Additional molecular properties.
  - Hartree-Fock and MP2 energy calculation.

Documentation for these functions is provided at the beginning of the `main.c` file.

## How to Set Up and Run

To install dependencies, compile, and run the program, follow the steps described in the [INSTALL.md](INSTALL.md) file.

### Quick Start

Once you have completed the setup, you can run the program using the script:
```bash
./compile_and_run.sh <file_path> <molecule>
```
Replace `<file_path>` with the path to the `.h5` input file and `<molecule>` with the molecule name. 

### Example Usage

The H2O molecule can be run as follows:

```bash
./compile_and_run.sh ./data/h2o.h5 H2O
```

### Testing

You can test the program with the `test.sh` script
```bash
./test.sh
```
This script will compile and run the program with the `.h5` provided input files.

