# Installation

## Prerequisites

The program depends on the following libraries:
- `HDF5`
- `TREXIO`

### Installing HDF5
On Ubuntu, you can install the `HDF5` library with:
```bash
sudo apt-get install libhdf5-dev
```

### Installing TREXIO
```bash
wget https://github.com/TREX-CoE/trexio/releases/download/v2.5.0/trexio-2.5.0.tar.gz
tar -zxvf trexio-2.5.0.tar.gz
cd trexio-2.5.0
./configure
make
sudo make install
```


### Compilation and running the program
Clone the repository and compile the program using the provided script:
```bash
git clone https://github.com/lavinia-santos/tccm-homeworks.git
cd project1/
./compile_and_run.sh <file_path> <molecule>
```

Replace `<file_path>` with the path to the `.h5` input file and `<molecule>` with the molecule name. The executable will be generated in the `lib` directory.



