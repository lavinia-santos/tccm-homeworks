/**
 * \file main.c
 * \brief This program reads the '.h5' input file and computes the HF energy along with MP2 correction.
 * 
 * \author Dijana Mitrovic, Lavínia Gabriela Teodoro dos Santos
 * \version 1.0
 * \date 2025-09-01
 * \copyright GNU Public License V3.0
 * This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    The authors can be contacted through email:
    - Dijana Mitrovic: dijana.mitrovic@univ-tlse3.fr
    - Lavínia Gabriela Teodoro dos Santos: lavinia-gabriela.teodoro-dos-santos@univ-tlse3.fr
 */


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
 * @param more_info [IN] flag to include additional results
 * @return double: the repulsion energy
 */
double compute_repulsion_energy(trexio_t *trexio_file, char *molecule, bool more_info);

/**
 * @brief Computes the number of occupied orbitals of the molecule associated to the given file
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @param verbose [IN] flag to enable detailed logging
 * @return int32_t: the number of occupied orbitals
 */
int32_t compute_occupied_orbitals(trexio_t *trexio_file, char *molecule, bool verbose);

/**
 * @brief Computes the number of molecular orbitals of the molecule associated to the given file
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @param verbose [IN] flag to enable detailed logging
 * @return int32_t: the number of molecular orbitals
 */
int32_t compute_number_of_molecular_orbitals(trexio_t *trexio_file, char *molecule, bool verbose);

/**
 * @brief Computes the one electron integrals of the molecule associated to the given file
 *
 * LLR:
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @param verbose [IN] flag to enable detailed logging
 * @return double*: the one electron integrals
 */
double *compute_one_electron_integrals(trexio_t *trexio_file, char *molecule, bool verbose);

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

/**
 * @brief Computes the flattened index for accessing two-electron integrals
 *
 * @param i [IN] index of the first molecular orbital
 * @param j [IN] index of the second molecular orbital
 * @param k [IN] index of the third molecular orbital
 * @param l [IN] index of the fourth molecular orbital
 * @param mo_num [IN] total number of molecular orbitals
 * @return int64_t: the flattened index for the given orbital combination
 */
int64_t get_index(int i, int j, int k, int l, int mo_num);

/**
 * @brief Computes the Hartree-Fock energy for a given molecule
 *
 * @param trexio_file [IN] the file where the data is stored
 * @param molecule [IN] the molecule considered
 * @param repulsion_energy [IN] nuclear repulsion energy of the molecule
 * @param one_e_integral [IN] one-electron integrals of the molecule
 * @param mo_num [IN] total number of molecular orbitals
 * @param more_info [IN] flag to include additional results
 * @param verbose [IN] flag to enable detailed logging
 * @return double: Hartree-Fock energy
 */
double compute_HF_energy(trexio_t *trexio_file, char *molecule, double repulsion_energy, double *one_e_integral, int mo_num, bool more_info, bool verbose);

/**
 * @brief Reads the molecular orbital energies from the given TREXIO file
 *
 * @param file [IN] the TREXIO file from which the data is read
 * @param mo_energy [OUT] array to store the molecular orbital energies
 * @return trexio_exit_code: status code indicating success or failure of the operation
 */
trexio_exit_code trexio_read_mo_energy(trexio_t *const file, double *const mo_energy);

/**
 * @brief Computes the MP2 energy correction for a given molecule
 *
 * @param buffer_size [IN] the number of read integrals
 * @param index [IN] the indexes of the two-electron integrals
 * @param value [IN] the values of the two-electron integrals
 * @param mo_num [IN] total number of molecular orbitals
 * @param mo_occ [IN] number of occupied molecular orbitals
 * @param orbital_energies [IN] orbital energy values
 * @return double: MP2 energy correction
 */
double compute_mp2_correction(int64_t buffer_size, int32_t *index, double *value, int mo_num, int mo_occ, double *orbital_energies);


/**
 * @brief Displays a short license notice.
 */
void display_interactive_license_notice();

/* PUBLIC FUNCTIONS DEFINITIONS */
void main(int argc, char *argv[])
{
    display_interactive_license_notice();
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
    // choose whether to print more information
    bool more_info = false;
    //choose whether to be verbose
    bool verbose = false;


    /* Compute Nuclear Repulsion Energy */
    double energy;
    energy = compute_repulsion_energy(trexio_file, molecule, more_info);

    /* Compute the number of occupied orbitals */
    int32_t n_up;
    n_up = compute_occupied_orbitals(trexio_file, molecule, verbose);

    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    mo_num = compute_number_of_molecular_orbitals(trexio_file, molecule, verbose);

    /* Compute the one electron integrals */
    double *data;
    data = compute_one_electron_integrals(trexio_file, molecule, verbose);

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
        if (i <= 10 && verbose)
        {
            printf("%dth integral, corresponding to <%d %d|%d %d> for molecule %s: %f\n", n, i, j, k, l, molecule, value[n]);
        }
    }
    // Compute HF energy
    double hf_energy = compute_HF_energy(trexio_file, molecule, energy, data, mo_num, more_info, verbose);

    /* Read molecular orbital energies */
    double *mo_energy = malloc(mo_num * sizeof(double));
    if (mo_energy == NULL)
    {
        printf("Memory allocation error for mo_energy\n");
        exit(1);
    }

    rc = trexio_read_mo_energy(trexio_file, mo_energy);
    if (rc != TREXIO_SUCCESS)
    {
        printf("TREXIO Error reading molecular orbital energies:\n%s\n",
               trexio_string_of_error(rc));
        free(mo_energy);
        exit(1);
    }
    else if (verbose)
    {
        for (int i = 0; i < mo_num; ++i)
        {
            printf("Molecular orbital energy %d: %.5f\n", i, mo_energy[i]);
        }
    }
    // Compute MP2 correction
    double mp2_correction = compute_mp2_correction(buffer_size, index, value, mo_num, n_up, mo_energy);
    printf("MP2 correction: %.8f\n", mp2_correction);

    // Compute total energy
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
double compute_repulsion_energy(trexio_t *trexio_file, char *molecule, bool more_info)
{
    trexio_exit_code rc;
    /* Compute Nuclear Repulsion Energy */
    double energy;
    rc = trexio_read_nucleus_repulsion(trexio_file, &energy);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        if (more_info)
        {
            printf("Nuclear repulsion energy of molecule %s: %.5f\n", molecule, energy);
        }
    }
    else
    {
        printf("TREXIO Error reading nuclear repulsion energy:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return energy;
}

int32_t compute_occupied_orbitals(trexio_t *trexio_file, char *molecule, bool verbose)
{
    trexio_exit_code rc;
    /* Compute the number of occupied orbitals */
    int32_t n_up;
    rc = trexio_read_electron_up_num(trexio_file, &n_up);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        if (verbose)
        {
            printf("Number of occupied orbitals of molecule %s: %d\n", molecule, n_up);
        }
    }
    else
    {
        printf("TREXIO Error reading number of occupied orbitals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return n_up;
}

int32_t compute_number_of_molecular_orbitals(trexio_t *trexio_file, char *molecule, bool verbose)
{
    trexio_exit_code rc;
    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    rc = trexio_read_mo_num(trexio_file, &mo_num);
    /* Check the return code to be sure reading was OK */
    if (rc == TREXIO_SUCCESS)
    {
        if (verbose)
        {
            printf("Number of molecular orbitals of molecule %s: %d\n", molecule, mo_num);
        }
    }
    else
    {
        printf("TREXIO Error reading number of molecular orbitals:\n%s\n",
               trexio_string_of_error(rc));
        exit(1);
    }
    return mo_num;
}

double *compute_one_electron_integrals(trexio_t *trexio_file, char *molecule, bool verbose)
{
    trexio_exit_code rc;

    /* Compute the number of molecular orbitals */
    int32_t mo_num;
    mo_num = compute_number_of_molecular_orbitals(trexio_file, molecule, verbose);

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
                if (i <= 10 && verbose)
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
    // free(data);
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


int64_t get_index(int i, int j, int k, int l, int mo_num)
{
// Calculate the index for accessing integrals
    return i * mo_num * mo_num * mo_num + j * mo_num * mo_num + k * mo_num + l;
}

double compute_HF_energy(trexio_t *trexio_file, char *molecule, double repulsion_energy, double *one_e_integral, int mo_num, bool more_info, bool verbose)
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
    // bool more_info = true;

    // Determine the number of occupied orbitals
    mo_occ = compute_occupied_orbitals(trexio_file, molecule, verbose);



    // Compute the sum of one-electron integrals for occupied orbitals
    for (i = 0; i < mo_occ; ++i)
    {
        sum_one_e += one_e_integral[i * mo_num + i];
    }
    if (more_info){
    printf("Sum of one-electron integrals for molecule %s: %.8f\nsum*2: %.8f\n", molecule, sum_one_e, sum_one_e * 2);
    }
    // Compute the sum of two-electron integrals for occupied orbitals
    compute_two_electrons_integrals(trexio_file, molecule, &buffer_size, &index, &value);
    // printf("Number of read two-electrons integrals for molecule %s: %ld\n", molecule, buffer_size);
    
    // Allocate memory for storing two-electron integral values
    size_t max_size = (size_t) mo_occ * mo_occ * mo_occ * mo_occ;
    double *two_e_value = (double *)calloc(max_size, sizeof(double));
    if (two_e_value == NULL)
    {
        printf("Error: Fail to allocate memory for two_e_value.\n");
        exit(1);
    }
     // Store two-electron integrals in a structured format
    for (n = 0; n < buffer_size; ++n)
    {
 
        i = index[4 * n + 0];
        j = index[4 * n + 1];
        k = index[4 * n + 2];
        l = index[4 * n + 3];
        // printf("i: %d, j: %d, k: %d, l: %d\n", i, j, k, l);

        if (i < mo_occ && j < mo_occ && k < mo_occ && l < mo_occ)
        {
            two_e_value[get_index(i, j, k, l, mo_occ)] = value[n];
            two_e_value[get_index(i, l, k, j, mo_occ)] = value[n];
            two_e_value[get_index(k, l, i, j, mo_occ)] = value[n];
            two_e_value[get_index(k, j, i, l, mo_occ)] = value[n];
            two_e_value[get_index(j, i, l, k, mo_occ)] = value[n];
            two_e_value[get_index(l, i, j, k, mo_occ)] = value[n];
            two_e_value[get_index(l, k, j, i, mo_occ)] = value[n];
            two_e_value[get_index(j, k, l, i, mo_occ)] = value[n];
        }
    }
    // printf("Two-electron integrals stored\n");
    // Compute energy contribution from two-electron integrals
    for (int i = 0; i < mo_occ; i++)
    {
        for (int j = 0; j < mo_occ; j++)
        {
            for (int k = 0; k < mo_occ; k++)
            {
                for (int l = 0; l < mo_occ; l++)
                {
                    if (i == k && j == l)
                    {
                        double ijij = two_e_value[get_index(i, j, i, j, mo_occ)];
                        double ijji = two_e_value[get_index(i, j, j, i, mo_occ)];
                        sum_two_e += 2 * ijij - ijji;
                        if (verbose){
                            printf("ijij i: %d, j: %d, k: %d, l: %d, value:%f \n", i, j, k, l, two_e_value[get_index(i, j, i, j, mo_num)]);
                            printf("ijji i: %d, j: %d, k: %d, l: %d, value:%f \n", i, j, j, i, two_e_value[get_index(i, j, j, i, mo_num)]);
                        }
                    }
                }
            }
        }
    }
    free(two_e_value);


    if (more_info){

    printf("Sum of two-electron integrals for molecule %s: %.7f\n", molecule, sum_two_e);
    }
    // Compute Hartree-Fock energy
    hf_energy = repulsion_energy + 2 * sum_one_e + sum_two_e;

    // Output the computed HF energy
    printf("Computed HF energy for molecule %s: %.7f\n", molecule, hf_energy);

    return hf_energy;
}

double compute_mp2_correction(int64_t buffer_size, int32_t *index, double *value, int mo_num, int mo_occ, double *orbital_energies)
{

    double mp2_correction = 0.0;

    // Allocate memory for storing two-electron integral values
    size_t max_size = (size_t) mo_num * mo_num * mo_num * mo_num;
    double *two_e_value = (double *)calloc(max_size, sizeof(double));
    if (two_e_value == NULL)
    {
        printf("Error: Fail to allocate memory for two_e_value.\n");
        exit(1);
    }

    // Store MP2-relevant two-electron integrals
    for (int n = 0; n < buffer_size; ++n)
    {
        int a = index[4 * n + 0];
        int b = index[4 * n + 1];
        int i = index[4 * n + 2];
        int j = index[4 * n + 3];

        // Filter for MP2 integrals < a b | i j >
        if (a >= mo_occ || b >= mo_occ || i < mo_occ || j < mo_occ)
        {
            two_e_value[get_index(a, b, i, j, mo_num)] = value[n];
            two_e_value[get_index(a, j, i, b, mo_num)] = value[n];
            two_e_value[get_index(i, j, a, b, mo_num)] = value[n];
            two_e_value[get_index(i, b, a, i, mo_num)] = value[n];
            two_e_value[get_index(b, a, j, i, mo_num)] = value[n];
            two_e_value[get_index(j, a, b, i, mo_num)] = value[n];
            two_e_value[get_index(j, i, b, a, mo_num)] = value[n];
            two_e_value[get_index(b, i, j, a, mo_num)] = value[n];
        }
        else
        {
            continue;
        }
    }

    // Compute MP2 energy correction
    for (int i = 0; i < mo_occ; i++)
    {
        for (int j = 0; j < mo_occ; j++)
        {
            for (int a = mo_occ; a < mo_num; a++)
            {
                for (int b = mo_occ; b < mo_num; b++)
                {
                    double e_i = orbital_energies[i];
                    double e_j = orbital_energies[j];
                    double e_a = orbital_energies[a];
                    double e_b = orbital_energies[b];
                    double ijab = two_e_value[get_index(a, b, i, j, mo_num)];
                    double ijba = two_e_value[get_index(b, a, i, j, mo_num)];
                    double denominator = e_i + e_j - e_a - e_b;
                    double numerator = ijab * (2 * ijab - ijba);
                    mp2_correction += numerator / denominator;
                }
            }
        }
    }
    free(two_e_value);
    return mp2_correction;
}



void display_interactive_license_notice() {
    printf("\n\nComputation of the MP2 energy version 1.0, Copyright (C) 2025\n");
    printf("This program comes with ABSOLUTELY NO WARRANTY.\n");
    printf("This is a free software, and you are welcome to redistribute it\n");
    printf("under certain conditions. For details, see the LICENSE file.\n\n\n");
}

