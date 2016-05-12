#include "Common.h"
#include "SearchAlgs.h"
#include "partition.h"

// Run multiple-chain Markov Chain Monte Carlo to partition nodes
// of a network into a set number of groups, with the option to
// correct for degree
//
// Author: Elizabeth Sander (esander@uchicago.edu)
// Version: 1.0.0
//
// Requirements: gcc compiler, GSL and GSL CBLAS libraries for C
//
// Inputs:
// - number of species
// - file path to an adjacency matrix
// - random (integer) seed
// - number of steps for each chain to take
// - number of chains
// - maximum number of groups
// - 0/1 flag for degree correction
// - alpha parameter for gamma prior for the rate parameter omega
// (degree corrected model only)
// - beta parameter for gamma prior for the rate parameter omega
// (degree corrected model only)
//
// Example function call from terminal, if NOT correcting for degree:
// ./FindGroups N ./path/to/file seed steps nchains maxgroups degreeflag
// Example function call from terminal, if correcting for degree:
// ./FindGroups N ./path/to/file seed steps nchains maxgroups degreeflag alpha beta
//
// A chain swap is attempted every 20 steps. This is hardcoded in the
// global variables in SearchAlgs.c
//
// Prints a vector of length N to file. This vector contains the best
// partition found by the algorithm. Each integer i in the vector
// represents the group identity of species i.
// The output file will be named:
// "<./path/to/file>-G-<maxgroups>-DC-<degreeflag>-alpha-<alpha>-beta-<beta>-Marginal<marginal loglikelihood>.txt"
//
// Therefore, you should copy your adjacency matrix to the folder where
// you want your output file.
// Note that if not using the degree-corrected model, <alpha> and
// <beta> in the output file name will be NA.

int main (int argc, char *argv[]){
    int N = atoi(argv[1]); // number of species
    char * FileName = argv[2]; // file storing the adjacency matrix
    // (if not correcting for degree, this matrix may be signed)
    int seed = atoi(argv[3]); // random seed
    int Steps = atoi(argv[4]); // number of steps
    int nChains = atoi(argv[5]); // for MC3, this is the number of
    // chains
    int maxGroups = atoi(argv[6]); // maximum number of groups allowed
    int DegreeCorrection = atoi(argv[7]); // 0 = no correction, 1 =
                                          // with correction
    double alpha, beta;
    if(DegreeCorrection == 1){
        alpha = atof(argv[8]); // for degree correction, alpha
                                      // parameter
        beta = atof(argv[9]); // for degree correction, beta
                                      // parameter
    }

    fprintf(stderr, "starting MC3 with %d steps, %d species, %d chains, %d groups, and degree correction %d\n", Steps, N, nChains, maxGroups, DegreeCorrection);

    // Read the adjacency matrix
    gsl_matrix_short * Adj = gsl_matrix_short_calloc(N,N);
    FILE * F;
    F=fopen(FileName, "rb");
    gsl_matrix_short_fscanf(F, Adj);
    fclose(F);

    // Set the random seed
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set (r, seed);

    // Best Solution
    double MarginalLikelihood = 0.0;
    // Stuff for output
    char OutFileName[1000];
    gsl_vector_short * BestSolution = gsl_vector_short_calloc(N);

    int maxLog = 2; //used for lookup table later
    int i,j;
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            if(gsl_matrix_short_get(Adj, i, j) != 0){
                maxLog++;
            }
        }
    }
    //set up log gamma function lookup table
    //max lgamma is 3 + N*N for the signed group model
    //for the degree corrected model, it's N + L
    int maxlgamma;
    if(DegreeCorrection == 0){
        maxlgamma = 3 + (N*N);
    }
    else {
        maxlgamma = N + maxLog - 2;
    }
    gsl_vector * lgammaLookup = gsl_vector_calloc(maxlgamma);
    for(i=0; i<maxlgamma; i++){
        gsl_vector_set(lgammaLookup, i, gsl_sf_lngamma(i+1));
    }
    gsl_vector * logLookup = gsl_vector_calloc(maxLog);
    for(i=0; i<maxLog; i++){
        gsl_vector_set(logLookup, i, log(i+1));
    }

    // Run the search
    MarginalLikelihood = MC3(N, Adj, Steps, nChains,
                             BestSolution, r, lgammaLookup, logLookup,
                             maxGroups, DegreeCorrection,
                             alpha, beta);
    fprintf(stderr, "MC3 complete. Best solution likelihood %.4f\n", MarginalLikelihood);
  
    // Save the results
    if(DegreeCorrection == 0){
        sprintf(OutFileName,"%s-G-%d-DC-%d-alpha-NA-beta-NA-Marginal%f", FileName,
                maxGroups, DegreeCorrection, MarginalLikelihood);
    }
    else {
        sprintf(OutFileName,"%s-G-%d-DC-%d-alpha-%f-beta-%f-Marginal%f",
                FileName, maxGroups, DegreeCorrection,
                alpha, beta, MarginalLikelihood);
    }
    F=fopen(OutFileName,"w");
    // print best solution
    PrintVectorShort(F, BestSolution);
    // close file
    fclose(F);
    // Free memory
    gsl_matrix_short_free(Adj);
    gsl_vector_free(lgammaLookup);
    gsl_vector_free(logLookup);
    gsl_vector_short_free(BestSolution);
    gsl_rng_free(r);
    return 0;
}
