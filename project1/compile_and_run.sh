if [ "$#" -ne 2 ]
then
    # Not enough input arguments
    echo "Usage compile_and_run.sh <file_path> <molecule>"
    echo "<file_path> The TREXIO file containing the data for the molecule"
    echo "<molecule> The molecule used for the computations (for print purposes)"
    echo "Exit"
    exit
fi

file_path=$1
molecule=$2

mkdir -p lib
gcc -I/usr/local/include src/main.c -L/usr/local/lib -ltrexio -o ./lib/main
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
./lib/main ${file_path} ${molecule}