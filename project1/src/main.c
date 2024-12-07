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
        printf("%dth integral, corresponding to <%d %d|%d %d> for molecule %s: %f\n", n, i, j, k, l, molecule, value[n]);
    }

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
                printf("One electron integrals of molecule %s: element (%d,%d) %.5f\n", molecule, i, j, data[i * mo_num + j]);
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