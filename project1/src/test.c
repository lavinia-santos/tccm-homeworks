#include <trexio.h>
#include <stdlib.h>
#include <stdio.h>

/* PRIVATE FUNCTIONS DECLARATIONS */
/**
 * @brief Computes the repulsion energy of the molecule associated to the given file
 * 
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return double: the repulsion energy 
 */
double compute_repulsion_energy(trexio_t *trexio_file, char* molecule);

/**
 * @brief Computes the number of occupied orbitals of the molecule associated to the given file
 * 
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return int32_t: the number of occupied orbitals
 */
int32_t compute_occupied_orbitals(trexio_t *trexio_file, char* molecule);

/**
 * @brief Computes the number of molecular orbitals of the molecule associated to the given file
 * 
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return int32_t: the number of molecular orbitals
 */
int32_t compute_number_of_molecular_orbitals(trexio_t *trexio_file, char* molecule);

/**
 * @brief Computes the one electron integrals of the molecule associated to the given file
 * 
 * LLR:
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return double*: the one electron integrals
 */
double* compute_one_electron_integrals(trexio_t *trexio_file, char* molecule);

/**
 * @brief Computes the two electrons integrals of the molecule associated to the given file
 * 
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @param buffer_size [OUT] the number of read integrals
 * @param index [OUT] the indexes
 * @param value [OUT] the values of the integrals
 */
void compute_two_electrons_integrals(trexio_t *trexio_file, char* molecule, int64_t* buffer_size, int32_t** index, double** value);

double compute_HF_energy(trexio_t *trexio_file, char* molecule, double repulsion_energy, double* one_e_integral, int mo_num);

/* PUBLIC FUNCTIONS DEFINITIONS */
void main(int argc, char *argv[])
{
    /* Read arguments */
    char *filePath = argv[1];
    char *molecule = argv[2];

    /* Iterators for printing */
    int i;
    int j;
    int k;
    int l;
    int n;

    /* Open the file */
    trexio_exit_code rc;
    trexio_t *trexio_file = trexio_open(filePath, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS)
    {
        printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    /* Compute Nuclear Repulsion Energy */
    double energy;
    energy = compute_repulsion_energy(trexio_file, molecule);

    /* Compute the number of occupied orbitals */
    int32_t n_up;
    n_up = compute_occupied_orbitals(trexio_file, molecule);

    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    mo_num = compute_number_of_molecular_orbitals(trexio_file, molecule);

    /* Compute the one electron integrals */
    double *data;
    data = compute_one_electron_integrals(trexio_file, molecule);

    /* Read two-electron integrals */
    int64_t buffer_size;
    int32_t* index;
    double* value;
    compute_two_electrons_integrals(trexio_file, molecule, &buffer_size, &index, &value);
    for (n=0; n<buffer_size; ++n)
    {
        i = index[4*n+0];
        j = index[4*n+1];
        k = index[4*n+2];
        l = index[4*n+3];
	if (i<=4)
	{
		printf("%dth integral, corresponding to <%d %d|%d %d> for molecule %s: %f\n", n, i, j, k, l, molecule, value[n]);
	}
    }
    double hf_energy;
    hf_energy = compute_HF_energy(trexio_file, molecule, energy, data, mo_num);
    printf("HF energy: %.5f\n", hf_energy);
    	
    /* Close the file */
    rc = trexio_close(trexio_file);
    if (rc != TREXIO_SUCCESS)
    {
        printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }
    trexio_file = NULL;
}

/* PRIVATE FUNCTIONS DEFINITIONS */
double compute_repulsion_energy(trexio_t *trexio_file, char* molecule)
{
    trexio_exit_code rc;
    /* Compute Nuclear Repulsion Energy */
    double energy;
    rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        printf("Nuclear repulsion energy of molecule %s: %.5f\n", molecule, energy);
    }
    else
    {
        printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return energy;
}

int32_t compute_occupied_orbitals(trexio_t *trexio_file, char* molecule)
{
    trexio_exit_code rc;
    /* Compute the number of occupied orbitals */
    int32_t n_up;
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        printf("Number of occupied orbitals of molecule %s: %d\n", molecule, n_up);
    }
    else
    {
        printf("TREXIO Error reading number of occupied orbitals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return n_up;
}

int32_t compute_number_of_molecular_orbitals(trexio_t *trexio_file, char* molecule)
{
    trexio_exit_code rc;
    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        printf("Number of molecular orbitals of molecule %s: %d\n", molecule, mo_num);
    }
    else
    {
        printf("TREXIO Error reading number of molecular orbitals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return mo_num;
}

double* compute_one_electron_integrals(trexio_t *trexio_file, char* molecule)
{
    trexio_exit_code rc;

    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    mo_num = compute_number_of_molecular_orbitals(trexio_file, molecule);

    /* Compute the one electron integrals */
    double* data;
    data = malloc(mo_num * mo_num * sizeof(double));
    rc = trexio_read_mo_1e_int_core_hamiltonian(trexio_file, data);

    /* iterators for printing */
    int i;
    int j;

    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        for (i = 0; i < mo_num; ++i)
        {
            for (j = 0; j < mo_num; ++j)
            {
		if (i <= 4)
		{
			printf("One electron integrals of molecule %s: element (%d,%d) %.5f\n", molecule, i, j, data[i * mo_num + j]);
		}
            }
        }
    }
    else
    {
        printf("TREXIO Error reading one electron integrals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return data;
}

void compute_two_electrons_integrals(trexio_t *trexio_file, char* molecule, int64_t* buffer_size, int32_t** index, double** value)
{
    trexio_exit_code rc;

    /* Read two-electron integrals */
    /* 1st step: Get the number of non-zero integrals */
    int64_t n_integrals;
    int64_t offset_file;
    rc = trexio_read_mo_2e_int_eri_size(trexio_file, &n_integrals);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        printf("Number of non-zero integrals for molecule %s: %ld\n", molecule, n_integrals);
        *index = malloc(4 * n_integrals * sizeof(int32_t));
        if (index == NULL)
        {
            fprintf(stderr, "Malloc failed for index");
            exit(1);
        }
        *value = malloc(n_integrals * sizeof(double));
        if (value == NULL)
        {
            fprintf(stderr, "Malloc failed for value");
            exit(1);
        }

        /* 2nd step: read the integrals from the file */
        offset_file = 0;
        *buffer_size = n_integrals;
        rc = trexio_read_mo_2e_int_eri(trexio_file, offset_file, buffer_size, *index, *value);
        /* Check the return code to be sure reading was OK */
        if (rc == TREXIO_SUCCESS)
        {
            printf("Number of read two electrons integrals for molecule %s: %ld\n", molecule, *buffer_size);
        }
        else
        {
            printf("TREXIO Error reading the two electrons integrals:\n%s\n",
               trexio_string_of_error(rc));
        }
    }
    else
    {
        printf("TREXIO Error reading the number of non-zero integrals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
}


struct permutation {
            int i;
            int j;
            int k;
            int l;
            float value;
        };

double compute_HF_energy(trexio_t *trexio_file, char* molecule, double repulsion_energy, double* one_e_integral, int mo_num) {
    trexio_exit_code rc;
    double sum_one_e = 0.0;  // Sum of one-electron integrals
    double sum_two_e = 0.0;  // Sum of two-electron integrals
    double hf_energy;
    int i, j, k, l, n;
    int mo_occ;  // Number of occupied orbitals
    int64_t buffer_size;
    int32_t *index;
    double *value;

    // Step 1: Determine the number of occupied orbitals
    mo_occ = compute_occupied_orbitals(trexio_file, molecule);
    printf("Number of occupied orbitals: %d\n", mo_occ);

    // Step 2: Compute the sum of one-electron integrals for occupied orbitals
    for (i = 0; i < mo_occ; ++i) {
        for (j = 0; j < mo_occ; ++j) {
            sum_one_e += one_e_integral[i * mo_num + j];
        }
    }

    // Step 3: Compute the sum of two-electron integrals for occupied orbitals
    compute_two_electrons_integrals(trexio_file, molecule, &buffer_size, &index, &value);
    printf("Number of read two electrons integrals for molecule %s: %ld\n", molecule, buffer_size);
    
    struct permutation perm[8*buffer_size];
    
    for (n = 0; n < buffer_size; ++n) {
        i = index[4 * n + 0];
        j = index[4 * n + 1];
        k = index[4 * n + 2];
        l = index[4 * n + 3];


        /*//create a matrix to write the permutations

        // int permutations[4][4];
            // if (i < mo_occ && j < mo_occ && k < mo_occ && l < mo_occ) {
                // printf("Permutations for <%d %d|%d %d>:\n", i, j, k, l);
                // printf("<%d %d|%d %d>: %d\n", i, l, k, j);
                // printf("<%d %d|%d %d>: %d\n", k, l, i, j);
                // struct permutation perm[4];

                // if ((k == i && l == j)) {
                // int count = 0;
                // write first combination ijkl to the matrix*/
            
                perm[8 * n+0].i = i;
                perm[8 * n+0].j = j;
                perm[8 * n+0].k = k;
                perm[8 * n+0].l = l;
                perm[8 * n+0].value = value[n];

                float v = value[n];
                if (i==j && j==k && k==l){
                    v = 0;
                }

                perm[8 * n+1].i = i;
                perm[8 * n+1].j = l;
                perm[8 * n+1].k = k;
                perm[8 * n+1].l = j;
                perm[8 * n+1].value = v;

                perm[8 * n+2].i = k;
                perm[8 * n+2].j = l;
                perm[8 * n+2].k = i;
                perm[8 * n+2].l = j;
                perm[8 * n+2].value = v;

                perm[8 * n+3].i = k;
                perm[8 * n+3].j = j;
                perm[8 * n+3].k = i;
                perm[8 * n+3].l = l;
                perm[8 * n+3].value = v;

                perm[8 * n+4].i = j;
                perm[8 * n+4].j = i;
                perm[8 * n+4].k = l;
                perm[8 * n+4].l = k;
                perm[8 * n+4].value = v;

                perm[8 * n+5].i = l;
                perm[8 * n+5].j = i;
                perm[8 * n+5].k = j;
                perm[8 * n+5].l = k;
                perm[8 * n+5].value = v;

                perm[8 * n+6].i = l;
                perm[8 * n+6].j = k;
                perm[8 * n+6].k = j;
                perm[8 * n+6].l = i;
                perm[8 * n+6].value = v;   

                perm[8 * n+7].i = j;
                perm[8 * n+7].j = k;
                perm[8 * n+7].k = l;
                perm[8 * n+7].l = i;
                perm[8 * n+7].value = v;          



                // printf("<%d %d|%d %d>: %f\n", perm[p].i, perm[p].j, perm[p].k, perm[p].l, value[n]);
                /*permutations[count][0] = i;
                permutations[count][1] = j;
                permutations[count][2] = k;
                permutations[count][3] = l;
                count++;
                // write the other combinations to the matrix
                permutations[count][0] = j;
                permutations[count][1] = i;
                permutations[count][2] = l;
                permutations[count][3] = k;
                count++;
                permutations[count][0] = l;
                permutations[count][1] = k;
                permutations[count][2] = j;
                permutations[count][3] = i;
                count++;
                permutations[count][0] = k;
                permutations[count][1] = l;
                permutations[count][2] = i;
                permutations[count][3] = j;*/
                for (int p = 0; p < 8; ++p) {
                    // perm[4 * n+p].i = i;
                    // perm[4 * n+p].j = j;
                    // perm[4 * n+p].k = k;
                    // perm[4 * n+p].l = l;
                    // perm[4 * n+p].value = value[n];

                    // printf("<%d %d|%d %d>: %f\n", perm[4*n+p].i, perm[4*n+p].j, perm[4*n+p].k, perm[4*n+p].l, perm[4*n+p].value);

                }      
                
                /*// create a matrix to write the permutations <ij|kl> as an int ijkl and the respective values

                // for (int p = 0; p < 4; ++p) {
                    
                
                //     printf("Term 2<ij|ij>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, k, l, value[n]);
                //     sum_two_e += 2 * value[n];
                //     //get permutation <ij|ji>
                //         printf("Permutation <ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, l, k, value[n]);
                //         sum_two_e -= value[n];
                    
                //     // printf("Permutation <ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, l, k, value[n]);
                //     // sum_two_e -= value[n];
                // } 
                // else if ((k == j && l == i) || (k == i && l == j)) {
                //     printf("Term -<ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, k, l, value[n]);
                //     sum_two_e -= value[n];
                // }*/
        //    }
        }
    
 

    

    for (int x = 0; x < 8* buffer_size; ++x)
    {
        if (perm[x].i < mo_occ && perm[x].j < mo_occ && perm[x].k < mo_occ && perm[x].l < mo_occ) {

        if (perm[x].i == perm[x].k && perm[x].j == perm[x].l)
        {
            sum_two_e += 2*perm[x].value;
            printf("Term 2<ij|ij>: %dth integral, <%d %d|%d %d>: %f\n", x, perm[x].i, perm[x].j, perm[x].k, perm[x].l, perm[x].value);
            if (perm[x].i == perm[x].j && perm[x].j == perm[x].l && perm[x].l == perm[x].k)           
            {
                sum_two_e -= perm[x].value;
                printf("Term -<ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", x, perm[x].i, perm[x].j, perm[x].k, perm[x].l, perm[x].value);

            } else {
                for (int y = 0; y < 8*buffer_size; y++){
                if (perm[y].i == perm[x].i && perm[y].j == perm[x].j && perm[y].k == perm[y].j && perm[y].l == perm[y].i){
                    sum_two_e -= perm[y].value;
                    printf("Term -<ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", x, perm[y].i, perm[y].j, perm[y].k, perm[y].l, perm[y].value);
                break;
                }
            }
            }
            
        }
        /*if (perm[x].i == perm[x].l && perm[x].j == perm[x].k)
        {
            sum_two_e -= perm[x].value;
            printf("Term -<ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", x, perm[x].i, perm[x].j, perm[x].k, perm[x].l, perm[x].value);
    }*/
    }
    }
    // if (i == k && j == l){
    //     printf("Term 2<ij|ij>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, k, l, value[n]);
    //     sum_two_e += 2 * value[n];
    //     //get permutation <ij|ji>
    //     printf("Permutation <ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, l, k, value[n]);
    //     sum_two_e -= value[n];
    // }

    // Step 4: Compute Hartree-Fock energy
    hf_energy = repulsion_energy + 2 * sum_one_e + sum_two_e;

    // Output the computed HF energy
    printf("Computed HF energy for molecule %s: %.5f\n", molecule, hf_energy);

    return hf_energy;
}
    


  /*  // if (i < mo_occ && j < mo_occ && k < mo_occ && l < mo_occ){
    // printf("Permutations for %dth integral: %d\n", n, count);
    // printf("Permutations for <%d %d|%d %d>:\n", i, j, k, l, count);
    // printf("<%d %d|%d %d>: %d\n", i, l, k, j, count);
    // printf("<%d %d|%d %d>: %d\n", k, l, i, j, count);
    // printf("<%d %d|%d %d>: %d\n", k, j, i, l, count);
    // printf("<%d %d|%d %d>: %d\n", j, i, l, k, count);
    // printf("<%d %d|%d %d>: %d\n", l, i, j, k, count);
    // printf("<%d %d|%d %d>: %d\n", l, k, j, i, count);
    // printf("<%d %d|%d %d>: %d\n", j, k, l, i, count);
    // }



//     for (int p = 0; p < count; ++p) {
//             i = permutations[p][0];
//             j = permutations[p][1];
//             k = permutations[p][2];
//             l = permutations[p][3];
//         if (i < mo_occ && j < mo_occ && k < mo_occ && l < mo_occ) {
//             if (k == i && l == j) {
//                 printf("Term 2<ij|ij>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, k, l, value[n]);
//                 sum_two_e += 2 * value[n];
//                 // create the exchange terms, because we know that ⟨ij|kl⟩ = ⟨il|kj⟩ = ⟨kl|ij⟩ = ⟨kj|il⟩ = ⟨ji|lk⟩ = ⟨li|jk⟩ = ⟨lk|ji⟩ = ⟨jk|li⟩
                

//             } else if (k == j && l == i) {
//                 // Caso <ij|ji>, subtrair uma vez
//                 printf("Term -<ij|ji>: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, k, l, value[n]);
//                 sum_two_e -= value[n];
//             } else {
//                 // Outras permutações que não aparecem diretamente no arquivo
//                 printf("Other permutation: %dth integral, <%d %d|%d %d>: %f\n", n, i, j, k, l, value[n]);
//             }
//         }
//     }
// }



//     // Step 4: Compute Hartree-Fock energy
//     hf_energy = repulsion_energy + 2*sum_one_e + sum_two_e;

//     // Step 5: Cleanup memory
//     // free(index);
//     // free(value);

//     // Output the computed HF energy
//     printf("Computed HF energy for molecule %s: %.5f\n", molecule, hf_energy);

//     return hf_energy;
// }*/


