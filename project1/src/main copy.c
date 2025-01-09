#include <trexio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* PRIVATE FUNCTIONS DECLARATIONS */
/**
 * @brief Computes the repulsion energy of the molecule associated to the given file
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return double: the repulsion energy
 */
double compute_repulsion_energy(trexio_t *trexio_file, char *molecule);

/**
 * @brief Computes the number of occupied orbitals of the molecule associated to the given file
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return int32_t: the number of occupied orbitals
 */
int32_t compute_occupied_orbitals(trexio_t *trexio_file, char *molecule);

/**
 * @brief Computes the number of molecular orbitals of the molecule associated to the given file
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return int32_t: the number of molecular orbitals
 */
int32_t compute_number_of_molecular_orbitals(trexio_t *trexio_file, char *molecule);

/**
 * @brief Computes the one electron integrals of the molecule associated to the given file
 *
 * LLR:
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @return double*: the one electron integrals
 */
double *compute_one_electron_integrals(trexio_t *trexio_file, char *molecule);

/**
 * @brief Computes the two electrons integrals of the molecule associated to the given file
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @param buffer_size [OUT] the number of read integrals
 * @param index [OUT] the indexes
 * @param value [OUT] the values of the integrals
 */
void compute_two_electrons_integrals(trexio_t *trexio_file, char *molecule, int64_t *buffer_size, int32_t **index, double **value);


double compute_HF_energy(trexio_t *trexio_file, char *molecule, double repulsion_energy, double *one_e_integral, int mo_num);

trexio_exit_code trexio_read_mo_energy(trexio_t* const file, double* const mo_energy);



double compute_MP2_energy(int64_t n_integrals, int32_t* index, double* value, int mo_num, int n_occ, double* orbital_energies);

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
    int32_t *index;
    double *value;
    compute_two_electrons_integrals(trexio_file, molecule, &buffer_size, &index, &value);
    for (n = 0; n < buffer_size; ++n)
    {
        i = index[4 * n + 0];
        j = index[4 * n + 1];
        k = index[4 * n + 2];
        l = index[4 * n + 3];
        if (i <= 10)
        {
            // printf("%dth integral, corresponding to <%d %d|%d %d> for molecule %s: %f\n", n, i, j, k, l, molecule, value[n]);
        }
    }
    // compute HF energy
    double hf_energy = compute_HF_energy(trexio_file, molecule, energy, data, mo_num);


    /* Read molecular orbital energies */
    double *mo_energy = malloc(mo_num * sizeof(double));
    if (mo_energy == NULL) {
        printf("Memory allocation error for mo_energy\n");
        exit(1);
    }

    rc = trexio_read_mo_energy(trexio_file, mo_energy);
    if (rc != TREXIO_SUCCESS) {
        printf("TREXIO Error reading molecular orbital energies:\n%s\n",
            trexio_string_of_error(rc));
        free(mo_energy);
        exit(1);
    } else {
        for (int i = 0; i < mo_num; ++i) {
            // printf("Molecular orbital energy %d: %.5f\n", i, mo_energy[i]);
        }
    }
    // Compute MP2 correction
    double mp2_correction = compute_MP2_energy(buffer_size, index, value, mo_num, n_up, mo_energy);
    printf("MP2 correction: %.8f\n", mp2_correction);

    double total_energy = hf_energy + mp2_correction;
    printf("Total energy for molecule %s: %.8f\n", molecule, total_energy);


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
double compute_repulsion_energy(trexio_t *trexio_file, char *molecule)
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

int32_t compute_occupied_orbitals(trexio_t *trexio_file, char *molecule)
{
    trexio_exit_code rc;
    /* Compute the number of occupied orbitals */
    int32_t n_up;
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        // printf("Number of occupied orbitals of molecule %s: %d\n", molecule, n_up);
    }
    else
    {
        printf("TREXIO Error reading number of occupied orbitals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return n_up;
}

int32_t compute_number_of_molecular_orbitals(trexio_t *trexio_file, char *molecule)
{
    trexio_exit_code rc;
    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        // printf("Number of molecular orbitals of molecule %s: %d\n", molecule, mo_num);
    }
    else
    {
        printf("TREXIO Error reading number of molecular orbitals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return mo_num;
}

double *compute_one_electron_integrals(trexio_t *trexio_file, char *molecule)
{
    trexio_exit_code rc;

    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    mo_num = compute_number_of_molecular_orbitals(trexio_file, molecule);

    /* Compute the one electron integrals */
    double *data;
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
                if (i <= 10)
                {
                    // printf("One electron integrals of molecule %s: element (%d,%d) %.5f\n", molecule, i, j, data[i * mo_num + j]);
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

void compute_two_electrons_integrals(trexio_t *trexio_file, char *molecule, int64_t *buffer_size, int32_t **index, double **value)
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
        // printf("Number of non-zero integrals for molecule %s: %ld\n", molecule, n_integrals);
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
            // printf("Number of read two electrons integrals for molecule %s: %ld\n", molecule, *buffer_size);
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

int get_index(int i, int j, int k, int l, int mo_num){

    return i*mo_num*mo_num*mo_num + j*mo_num*mo_num + k*mo_num + l;
}


double compute_HF_energy(trexio_t *trexio_file, char *molecule, double repulsion_energy, double *one_e_integral, int mo_num)
{
    trexio_exit_code rc;
    double sum_one_e = 0.0; // Sum of one-electron integrals
    double sum_two_e = 0.0; // Sum of two-electron integrals
    double hf_energy;
    int i, j, k, l, n;
    int mo_occ; // Number of occupied orbitals
    int64_t buffer_size;
    int32_t *index;
    double *value;

    // Determine the number of occupied orbitals
    mo_occ = compute_occupied_orbitals(trexio_file, molecule);
    // printf("Number of occupied orbitals: %d\n", mo_occ);
    size_t max_size = (size_t) mo_occ * mo_occ * mo_occ * mo_occ;
    double* mo_tei = (double*) calloc(max_size, sizeof(double));
    // Compute the sum of one-electron integrals for occupied orbitals
    for (i = 0; i < mo_occ; ++i)
    {
        sum_one_e += one_e_integral[i * mo_num + i];
    }
    printf("Sum of one-electron integrals for molecule %s: %.8f\n sum*2: %.8f\n", molecule, sum_one_e, sum_one_e*2);

    // Compute the sum of two-electron integrals for occupied orbitals
    compute_two_electrons_integrals(trexio_file, molecule, &buffer_size, &index, &value);
    // printf("Number of read two-electrons integrals for molecule %s: %ld\n", molecule, buffer_size);

    float energy = 0;

    /*// Store the two-electron integrals in the 1D array
    // for (n = 0; n < buffer_size; ++n)
    // {
    //     i = index[4 * n + 0];
    //     j = index[4 * n + 1];
    //     k = index[4 * n + 2];
    //     l = index[4 * n + 3];

    //     V(i, j, k, l) = value[n];
    //     printf("v[%d][%d][%d][%d]: %f\n", i, j, k, l, V(i, j, k, l));*/
    // }

    // Compute the energy contribution from two-electron integrals
    for (n = 0; n < buffer_size; ++n)
    {
        i = index[4 * n + 0];
        j = index[4 * n + 1];
        k = index[4 * n + 2];
        l = index[4 * n + 3];
        // printf("i: %d, j: %d, k: %d, l: %d\n", i, j, k, l);

        if (mo_tei == NULL) {
            printf("Error: Unable to allocate memory for mo_tei.\n");
            exit(1);
        }

        if (i < mo_occ && j < mo_occ && k < mo_occ && l < mo_occ)
        {
            mo_tei[get_index(i, j, k, l, mo_num)] = value[n];
            mo_tei[get_index(i, l, k, j, mo_num)] = value[n];
            mo_tei[get_index(k, l, i, j, mo_num)] = value[n];
            mo_tei[get_index(k, j, i, l, mo_num)] = value[n];
            mo_tei[get_index(j, i, l, k, mo_num)] = value[n];
            mo_tei[get_index(l, i, j, k, mo_num)] = value[n];
            mo_tei[get_index(l, k, j, i, mo_num)] = value[n];
            mo_tei[get_index(j, k, l, i, mo_num)] = value[n];

            /*if (i == j && j == k && k == l)
            {
                // energy += 2 * value[n] - value[n];
                energy += value[n];
            }
            else if (i == k && j == l && i != j)
            {
                // energy += 4 * (2 * value[n] - value[n]);
                energy += 4 * value[n];
            }
            else if ((i == j && l == k) || (i == l && j == k))
            {
                energy -= 2 * value[n];
            }*/
        }

    }
        for (int i = 0; i < mo_occ; i++) { // occupied orbitals
        for (int j = 0; j < mo_occ; j++) { // occupied orbitals
            for (int k = 0; k < mo_occ; k++) { // virtual orbitals
                for (int l = 0; l < mo_occ; l++) { // virtual orbitals
                if (i == k && j == l) {
                    //
                //  printf("ijij i: %d, j: %d, k: %d, l: %d, value:%f \n", i, j, k, l, mo_tei[get_index(i, j, i, j, mo_num)]);
                //  printf("ijji i: %d, j: %d, k: %d, l: %d, value:%f \n", i, j, j, i, mo_tei[get_index(i, j, j, i, mo_num)]);
                    double ijij = mo_tei[get_index(i, j, i, j, mo_num)];
                    double ijji = mo_tei[get_index(i, j, j, i, mo_num)];
                    sum_two_e += 2*ijij - ijji;
                }
                }
            }
        }
    }

    printf("Sum of two-electron integrals for molecule %s: %.8f\n", molecule, sum_two_e);
    


    // Step 4: Compute Hartree-Fock energy
    hf_energy = repulsion_energy + 2 * sum_one_e + sum_two_e;

    // Output the computed HF energy
    printf("Computed HF energy for molecule %s: %.7f\n", molecule, hf_energy);

    return hf_energy;
}




double compute_MP2_energy(int64_t buffer_size, int32_t* index, double* value, int mo_num, int n_occ, double* orbital_energies) {

    double MP2_energy = 0.0;
    

    size_t max_size = (size_t) mo_num * mo_num * mo_num * mo_num;
    double* mo_tei = (double*) calloc(max_size, sizeof(double));
    if (mo_tei == NULL) {
        printf("Error: Unable to allocate memory for mo_tei.\n");
        exit(1);
    }
    

    for (int n = 0; n < buffer_size; ++n) {
        int i = index[4*n + 0];
        int j = index[4*n + 1];
        int k = index[4*n + 2];
        int l = index[4*n + 3];

        // Filter only MP2 relevant integrals < a b | i j >
        if (i < n_occ || j < n_occ || k >= n_occ || l >= n_occ) {
        continue; 
        }
        mo_tei[get_index(i, j, k, l, mo_num)] = value[n];
        mo_tei[get_index(i, l, k, j, mo_num)] = value[n];
        mo_tei[get_index(k, l, i, j, mo_num)] = value[n];
        mo_tei[get_index(k, j, i, l, mo_num)] = value[n];
        mo_tei[get_index(j, i, l, k, mo_num)] = value[n];
        mo_tei[get_index(l, i, j, k, mo_num)] = value[n];
        mo_tei[get_index(l, k, j, i, mo_num)] = value[n];
        mo_tei[get_index(j, k, l, i, mo_num)] = value[n];

    }

    for (int i = 0; i < n_occ; i++) { // occupied orbitals
        for (int j = 0; j < n_occ; j++) { // occupied orbitals
            for (int a = n_occ; a < mo_num; a++) { // virtual orbitals
                for (int b = n_occ; b < mo_num; b++) { // virtual orbitals
                    double coulomb = mo_tei[get_index(a, b, i, j, mo_num)];
                    double exchange = mo_tei[get_index(b, a, i, j, mo_num)];
                    double denominator = orbital_energies[i] + orbital_energies[j] - orbital_energies[a] - orbital_energies[b];
                    double numerator = coulomb * (2*coulomb - exchange);
                    MP2_energy += numerator / denominator;
                }
            }
        }
    }
    free(mo_tei);
    return MP2_energy;
}








