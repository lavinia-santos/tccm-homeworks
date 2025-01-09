declare -a arr=("c2h2" "ch4" "co2" "h2o" "h3coh" "hcn")

# Loop through the array
for molecule in "${arr[@]}"
do
   ./compile_and_run.sh "data/${molecule}.h5" "$molecule"
done

