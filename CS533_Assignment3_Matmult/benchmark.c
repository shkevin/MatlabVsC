#include <stdlib.h>     // For: exit, drand48, malloc, free, NULL, EXIT_FAILURE
#include <stdio.h>      // For: perror
#include <string.h>     // For: memset

#include <float.h>      // For: DBL_EPSILON
#include <math.h>       // For: fabs

#include <sys/time.h>   // For struct timeval, gettimeofday
#include <sys/types.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <string.h>

#ifdef __INTEL_COMPILER
#include <mkl.h>
#else
#include <gsl/gsl_cblas.h>
#endif

#include "dgemm.h"
#include <time.h>       // For: time

extern char *optarg;
extern int optind;
extern int optopt;
extern int opterr;
extern int optreset;

#define MAX_SPEED (2.67*4)

/* reference_dgemm wraps a call to the BLAS-3 routine DGEMM, via the C
 * interface - hence the reference semantics. */
void reference_dgemm(int N, double ALPHA, double A[N][N], double B[N][N],
             double C[N][N])
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, ALPHA,
            (double *)A, N, (double *)B, N, 1.0, (double *)C, N);
}

void usage(char *name)
{
    fprintf(stderr, "usage: %s [-h] [-s size] [-o outfile]\n", name);
}

double wall_time()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return 1. * t.tv_sec + 1.e-6 * t.tv_usec;
}

void die(const char *message)
{
    perror(message);
    exit(EXIT_FAILURE);
}

void fill(double *p, int n)
{
    for (int i = 0; i < n; ++i)
        p[i] = 2 * drand48() - 1;   // Uniformly distributed over [-1, 1]
}

void absolute_value(double *p, int n)
{
    for (int i = 0; i < n; ++i)
        p[i] = fabs(p[i]);
}

int *process_sizes(char *arg, int *nitems, int *max)
{
    int num_items, max_val = -1, i;
    int *items = NULL;
    char *p, *curr;

    if (!arg)
        return NULL;

    /* First, count the number of sizes */
    num_items = 0;
    p = arg;
    while (*p) {
        while (*p && *p != ',')
            p++;
        num_items++;
        if (*p)
            p++;
    }

    /* If no items, bail. */
    if (!num_items)
        goto error;

    /* Now create an array of these items. */
    items = malloc(sizeof(int) * num_items);
    curr = p = arg;
    i = 0;
    while (*curr) {
        int val;
        /* Find the end of this item */
        val = strtol(curr, &curr, 10);
        if (val <= 0)
            goto error;
        if (val > max_val)
            max_val = val;
        items[i++] = val;
        if (*curr)
            curr++;
    }

    if (nitems)
        *nitems = num_items;
    if (max)
        *max = max_val;
    return items;

 error:
    if (items)
        free(items);
    *nitems = 0;
    *max = 0;
    return NULL;
}

int all_sizes[] =
 {31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
    319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769,
        1538, 1539, 3078, 3079, 6158, 6159, 6160, 10000};

// int default_sizes[] =
//  {31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
//     319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769,
//         1538, 1539, 3078, 3079, 6158, 6159, 6160, 10000};

int default_sizes[] =
 {31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
    319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767};

//Below code is from https://www.geeksforgeeks.org/shuffle-a-given-array-using-fisher-yates-shuffle-algorithm/
/**********************************************************************************/
// A utility function to print an array 
void printArray (int arr[], int n) 
{ 
    for (int i = 0; i < n; i++) 
        printf("%d ", arr[i]); 
    printf("\n"); 
} 

// A utility function to swap to integers 
void swap (int *a, int *b) 
{ 
    int temp = *a; 
    *a = *b; 
    *b = temp; 
} 
  
// A function to generate a random permutation of arr[] 
void randomize(int arr[], int n) 
{ 
    // Use a different seed value so that we don't get same 
    // result each time we run this program 
    srand (time(NULL)); 
  
    // Start from the last element and swap one by one. We don't 
    // need to run for the first element that's why i > 0 
    for (int i = n-1; i > 0; i--) 
    { 
        // Pick a random index from 0 to i 
        int j = rand() % (i+1); 
  
        // Swap arr[i] with the element at random index 
        swap(&arr[i], &arr[j]); 
    } 
} 
/**********************************************************************************/

int FindIndex(const int a[], size_t size, int value)
{
    int index = 0;

    while (index < size && a[index] != value ) ++index;

    return ( index == size ? -1 : index);
}

void create_folder(char *dirname)
{
    int check;
    check = mkdir(dirname,0777);

    if (check)
    {
      printf("Unable to create folder\n");
        exit(1);
    }
}

/* The benchmarking program */
int main(int argc, char **argv)
{
    int *test_sizes = NULL, nsizes, nmax;
    int numExperiment = 3;
    double cputime = 0;
    double walltime = 0;
    char cpu_filename[150];
    char wall_filename[150];

    test_sizes = default_sizes;
    nsizes = sizeof(default_sizes) / sizeof(default_sizes[0]);
    nmax = test_sizes[nsizes - 1];

    int iterations[] = {500,308,307,307,153,153,76,76,38,38,38,32,31,25,25,24,23,20,16,16,15,12,
        12,12,11,9,9,6,6,6,4,4,1,1};

    int iteration = 0;

    char parent[50];
    sprintf(parent, "./%s", argv[1]);
    create_folder(parent);

    char exp_str[75];
    sprintf(exp_str, "%s/%s", parent, "experimentInfo");
    create_folder(exp_str);

    char cpu_time[150];
    sprintf(cpu_time, "%s/%s", exp_str, "cputime");
    create_folder(cpu_time);

    char walltimeFold[150];
    sprintf(walltimeFold, "%s/%s", exp_str, "walltime");
    create_folder(walltimeFold);

    FILE *cpufile = NULL;
    FILE *wallfile = NULL;

    clock_t cpu_start, cpu_end;

    /* allocate memory for all problems */
    double *buf = NULL;
    buf = (double *)malloc(3 * nmax * nmax * sizeof(double));
    if (buf == NULL)
        die("failed to allocate largest problem size");


    //Randomize size locations
    randomize(test_sizes, nsizes);
    // printf("randomized array...\n");
    // for (int i = 0; i < nsizes; ++i)
    // {
    //     printf("%d, ", test_sizes[i]);
    // }
    // printf("\n");

    /* For each test size */
    for (int isize = 0; isize < nsizes; ++isize) {

        /* Create and fill 3 random matrices A,B,C */
        int n = test_sizes[isize];
        // printf("size: %d\n", n);

        double *A = buf + 0;
        double *B = A + nmax * nmax;
        double *C = B + nmax * nmax;

        fill(A, n * n);
        fill(B, n * n);
        fill(C, n * n);

        sprintf(cpu_filename, "%s/%d.txt", cpu_time,n);
        sprintf(wall_filename, "%s/%d.txt", walltimeFold,n);

        // printf("\n");
        cpufile = fopen(cpu_filename, "a");
        wallfile = fopen(wall_filename, "a");

        /* Warm-up cache for this size*/
        square_dgemm(n, (double (*)[])A, (double (*)[])B,
                 (double (*)[])C);

        // printf("%s\n", filename);

        /* Benchmark n_iterations runs of square_dgemm */
        iteration = iterations[FindIndex(all_sizes,nsizes,n)];
        // printf("iteration index = %d\n", iteration+1);
        // printf("iteration number = %d\n",iterations[iteration]);
        for (int i = 0; i < iteration; ++i)
        {
            for (int j = 0; j < numExperiment; ++j)
            {

                walltime = -wall_time();;
                cpu_start = clock();

                square_dgemm(n, (double (*)[])A,
                     (double (*)[])B, (double (*)[])C);

                cpu_end = clock();
                cputime = ((double) (cpu_end - cpu_start)) / CLOCKS_PER_SEC;
                walltime += wall_time();
                // printf("%lf    ", cputime);
                fprintf(cpufile, "%lf,",cputime);
                fprintf(wallfile, "%lf,",walltime);
            }
            // printf("\n");
            fprintf(cpufile,"\n");
            fprintf(wallfile,"\n");
        }
        cputime = 0;
        walltime = 0;
        fclose(cpufile);
        fclose(wallfile);
    }
    free(buf);

    return 0;
}
