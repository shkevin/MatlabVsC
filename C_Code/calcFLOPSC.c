#include <stdlib.h>     // For: exit, drand48, malloc, free, NULL, EXIT_FAILURE
#include <stdio.h>      // For: perror
#include <string.h>     // For: memset

#include <float.h>      // For: DBL_EPSILON
#include <math.h>       // For: fabs

#include <sys/time.h>   // For struct timeval, gettimeofday
#include <unistd.h>
#include <stdbool.h>    // For: bool
#include <time.h>       // For: time


int default_sizes[] =
 { 31, 32, 96, 97, 127, 128, 129, 191, 192, 229, 255, 256, 257,
    319, 320, 321, 417, 479, 480, 511, 512, 639, 640, 767, 768, 769};

/* ************************************************
* Purpose: Stores the commands given to program 
************************************************* */
typedef struct commands
{
    int opt, numCars, maxPerCar;
    bool verbose;
}commands;

void usage(char *name)
{
    fprintf(stderr, "usage: %s [-h] [-s size] [-o outfile]\n", name);
}

/* ************************************************
* PARAMETERS: None.
* FUNCTION: Prints the usage when the -h flag is
            set.
* RETURNS: None.
************************************************* */
void printUsage()
{
    printf("Usage: ./bin/PA04-KECO [-h] -v -N <num> -M <num>\n");
    printf("Options:\n");
    printf("  -h         Print this help message.\n");
    printf("  -v         Optional verbose flag.\n");
    printf("Examples:\n  linux> ./calcFLOPSC\n");
}

// A utility function to swap to integers 
void swap (int *a, int *b) 
{ 
    int temp = *a; 
    *a = *b; 
    *b = temp; 
} 

//Below code is from https://www.geeksforgeeks.org/shuffle-a-given-array-using-fisher-yates-shuffle-algorithm/
/**********************************************************************************/
// A utility function to print an array 
void printArray (int arr[], int n) 
{ 
    for (int i = 0; i < n; i++) 
        printf("%d ", arr[i]); 
    printf("\n"); 
} 
  
// A function to generate a random permutation of arr[] 
void randomize ( int arr[], int n ) 
{ 
    // Use a different seed value so that we don't get same 
    // result each time we run this program 
    srand ( time(NULL) ); 
  
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

commands parseCMD(int argc, char *argv[])
{
    commands com;

    while ((com.opt = getopt(argc, argv, "h:vN:M:")) != -1)
    {
        switch (com.opt)
        {
            case 'h':
            printUsage();
            break;
            case 'v':
            com.verbose = true;
            break;
            default:
            printf("Incorrect use\n");
            break;
        }
    }

    return com;
}

void fill(double *p, int n)
{
    for (int i = 0; i < n; ++i)
        p[i] = 2 * drand48() - 1;   // Uniformly distributed over [-1, 1]
}

int main(int argc, char **argv)
{
    FILE *outfile = NULL;
    commands com = parseCMD(argc,argv);

    int iterations = 10;
    int numExperiments = 10;
    int *sizes = &(default_sizes[0]);
    int len = sizeof(default_sizes)/sizeof(int);

    double *buf = NULL;
    int nmax = sizes[len - 1];
    buf = (double *)malloc(3 * nmax * nmax * sizeof(double));

    //Shuffle sizes
    randomize(sizes,len);

    for (int s = 0; s < len; ++s)
    {
        int n = sizes[s];
        int times[iterations];
        int iterationTime = 0;

        //Set times to all zeros
        memset(times, 0, sizeof times);
        int seconds = 0;

        /* Create and fill 3 random matrices A,B,C */
        double *A = buf + 0;
        double *B = A + nmax * nmax;
        double *C = B + nmax * nmax;

        fill(A, n * n);
        fill(B, n * n);
        fill(C, n * n);

        for (int i = 0; i < iterations; ++i)
        {
            

            for (int j = 0; j < numExperiments; ++j)
            {
                //Perform multiplication
            }


        }


    }

    free(buf);


    return 0;
}