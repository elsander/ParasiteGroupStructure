#include "Common.h"

//globals
int SWAPS=0;
int ACCEPTEDSWAPS=0;

//Printing functions-- mostly used for testing
int PrintMatrixShort(FILE * F, gsl_matrix_short * A){
    int NR = A->size1;
    int NC = A->size2;
    int i,j;
    for (i = 0;i < NR; i++){
        for (j = 0;j < NC; j++){
            fprintf(F, "%d ",gsl_matrix_short_get(A, i, j));
        }
        fprintf(F, "\n");
    }
    fprintf(F, "\n");
    return 0;
}

int PrintVectorFloat(FILE * F, gsl_vector * V){
    int N = V->size;
    int i;
    for (i = 0; i < N; i++){
        fprintf(F, "%1.4f ",gsl_vector_get(V, i));
    }
    fprintf(F, "\n");
    return 0;
}

int PrintVectorShort(FILE * F, gsl_vector_short * V){
    int N = V->size;
    int i;
    for (i = 0; i < N; i++){
        fprintf(F, "%d ",gsl_vector_short_get(V, i));
    }
    fprintf(F, "\n");
    return 0;
}

int Partition_Initialize(gsl_vector_short * V, int N, int maxGroups, gsl_rng * r){
    //Randomly generate the first grouping
    int i;
    //initialize randomly
    for (i = 0; i < N; i++){
        gsl_vector_short_set(V, i, gsl_rng_uniform_int(r, maxGroups));
    }

    return 0;
}

int RGF(int N, gsl_vector_short * Chain, gsl_vector_short * RGFswap){
    //RGF reorders groupings so that new groups are introduced in
    //order in the vector from 0 to N-1
    //Example Input vector: 2 3 0 2
    //Output vector: 0 1 2 0

    //clear RGFswap
    gsl_vector_short_set_all(RGFswap, -1);

    int i;
    int groupNum = 0;
    //RGF the entire vector as a whole unit
    for(i=0; i<N; i++){
        if(gsl_vector_short_get(RGFswap, gsl_vector_short_get(Chain, i)) == -1){
            //if there's a -1, we haven't seen that value yet
            gsl_vector_short_set(RGFswap, gsl_vector_short_get(Chain, i), groupNum);
            groupNum++;
        }
        //set the value from RGFswap to the ith element of Chain
        gsl_vector_short_set(Chain, i,
                             gsl_vector_short_get(RGFswap,
                                                  gsl_vector_short_get(Chain, i)));
    }
    return 0;
}

double Marginal(gsl_vector_short * V,
                gsl_vector_short * VCopy,
                gsl_matrix_short * Adj,
                int N,
                gsl_vector * lgammaLookup,
                gsl_vector * logLookup){
    //Calculate marginal log likelihood of a partition
    int i, j;
    int Groupi, Groupj;
    int NumPartitions;

    gsl_vector_short_memcpy(VCopy, V);
    NumPartitions = gsl_vector_short_max(VCopy) + 1;
    // store the information on the links between two groups
    int NonZeros[NumPartitions][NumPartitions];
    int Zeros[NumPartitions][NumPartitions];
    int Positives[NumPartitions][NumPartitions];
    int Negatives[NumPartitions][NumPartitions];

    // initialize
    for (i = 0; i < NumPartitions; i++){
        for (j = 0; j < NumPartitions; j++){
            NonZeros[i][j] = 0;
            Zeros[i][j] = 0;
            Positives[i][j] = 0;
            Negatives[i][j] = 0;
        }
    }
    // Assign link information
    for (i = 0; i < N; i++){
        Groupi = gsl_vector_short_get(VCopy, i);
        for (j = 0; j < N; j++){
            Groupj = gsl_vector_short_get(VCopy, j);
            if (gsl_matrix_short_get(Adj, i, j) != 0){
                NonZeros[Groupi][Groupj]++;
                if (gsl_matrix_short_get(Adj, i, j) > 0){
                    Positives[Groupi][Groupj]++;
                }
                else{
                    Negatives[Groupi][Groupj]++;
                }
            }
            else{
                Zeros[Groupi][Groupj]++;
            }
        }
    }

    // Compute Marginal Likelihood
    double tmp;
    double MarginalLikelihood = 0.0;
    for (i = 0; i < NumPartitions; i++){
        for (j = 0; j < NumPartitions; j++){
            //have to subtract one from index for lookups
            //since the 0th entry is lgamma(1), not lgamma(0)
            tmp = (gsl_vector_get(lgammaLookup, Positives[i][j]) +
                   gsl_vector_get(lgammaLookup, Zeros[i][j]) +
                   gsl_vector_get(lgammaLookup, Negatives[i][j]) -
                   gsl_vector_get(logLookup, NonZeros[i][j]) -
                   gsl_vector_get(lgammaLookup, 1 + NonZeros[i][j] + Zeros[i][j]));
            MarginalLikelihood += tmp;
        }
    }

    return (MarginalLikelihood);
}

double DegreeMarginal(gsl_vector_short * V,
                      gsl_vector_short * VCopy,
                      gsl_matrix_short * Adj,
                      int N,
                      double alpha,
                      double beta,
                      gsl_vector * lgammaLookup){
    //Calculate marginal log likelihood of a partition
    int i, j;
    int Groupi, Groupj;
    int NumPartitions;

    gsl_vector_short_memcpy(VCopy, V);
    NumPartitions = gsl_vector_short_max(VCopy) + 1;
    // store the information on the links between two groups
    int NonZeros[NumPartitions][NumPartitions];
    int Zeros[NumPartitions][NumPartitions];
    // number of species in each group
    int nr[NumPartitions];
    
    int indegree[N];
    int outdegree[N];
    int indegree_by_group[NumPartitions];
    int outdegree_by_group[NumPartitions];

    // initialize
    for(i = 0; i < N; i++){
        indegree[i] = 0;
        outdegree[i] = 0;
        for(j = 0; j < N; j++){
            indegree[i] += gsl_matrix_short_get(Adj, j, i);
            outdegree[i] += gsl_matrix_short_get(Adj, i, j);
        }
    }
    
    for (i = 0; i < NumPartitions; i++){
        for (j = 0; j < NumPartitions; j++){
            NonZeros[i][j] = 0;
            Zeros[i][j] = 0;
        }
        nr[i] = 0;
        indegree_by_group[i] = 0;
        outdegree_by_group[i] = 0;
    }
    // Assign link information
    for (i = 0; i < N; i++){
        Groupi = gsl_vector_short_get(VCopy, i);
        for (j = 0; j < N; j++){
            Groupj = gsl_vector_short_get(VCopy, j);
            if (gsl_matrix_short_get(Adj, i, j) != 0){
                NonZeros[Groupi][Groupj]++;
            }
            else{
                Zeros[Groupi][Groupj]++;
            }
        }
        nr[Groupi]++;
        indegree_by_group[Groupi] += indegree[i];
        outdegree_by_group[Groupi] += outdegree[i];
    }

    // Compute Marginal Likelihood
    double tmp;
    double MarginalLikelihood;
    // outer term
    MarginalLikelihood = NumPartitions * NumPartitions *
        (alpha*log(beta) - gsl_sf_lngamma(alpha));
    // sum over all species
    for(i = 0; i < N; i++){
        // equivalent to (indegree_i)!, because lgammaLookup
        // is shifted by 1 (lgamma[1] = gamma(2)),
        // and gamma(n) = (n-1)!
        MarginalLikelihood += gsl_vector_get(lgammaLookup, indegree[i]) +
            gsl_vector_get(lgammaLookup, outdegree[i]);
    }
    // sum over all blocks
    for (i = 0; i < NumPartitions; i++){
        for (j = 0; j < NumPartitions; j++){
            tmp = gsl_sf_lngamma(alpha + (double)NonZeros[i][j]) -
                (alpha - (double)NonZeros[i][j])*log(1 + beta);
            MarginalLikelihood += tmp;
        }
        //have to subtract one from index for lookups
        //since the 0th entry is lgamma(1), not lgamma(0)
        tmp = -2 * gsl_vector_get(lgammaLookup, nr[i] - 1) -
            gsl_vector_get(lgammaLookup, nr[i] + indegree_by_group[i] - 1) -
            gsl_vector_get(lgammaLookup, nr[i] + outdegree_by_group[i] - 1);
        MarginalLikelihood += tmp;
    }

    return (MarginalLikelihood);
}

double Partition_Likelihood(gsl_vector_short * V,
                            gsl_vector_short * VCopy,
                            gsl_matrix_short * Adj,
                            int N,
                            gsl_vector * lgammaLookup,
                            gsl_vector * logLookup,
                            int DegreeCorrection,
                            double alpha,
                            double beta){
    // call Bayes Factor or NML likelihood functions
    double likelihood;
    if(DegreeCorrection == 0){
        likelihood = Marginal(V, VCopy, Adj, N, lgammaLookup, logLookup);
    }
    else {
        likelihood = DegreeMarginal(V, VCopy, Adj, N, alpha, beta, lgammaLookup);
    }
    return(likelihood);
}

int TrySwap(int N,
            gsl_matrix_short * Adj,
            gsl_vector_short * HotChain,
            gsl_vector_short * ColdChain,
            gsl_vector_short * ChainCopy,
            double HotTemp, double ColdTemp,
            gsl_rng * r,
            gsl_vector * lgammaLookup,
            gsl_vector * logLookup,
            int DegreeCorrection,
            double alpha,
            double beta){

    double HotLik, ColdLik;
    double ProbSwap;

    //Attempt a swap between two chains

    HotLik = Partition_Likelihood(HotChain, ChainCopy,
                                  Adj, N, lgammaLookup, logLookup,
                                  DegreeCorrection,
                                  alpha, beta);
    ColdLik = Partition_Likelihood(ColdChain, ChainCopy,
                                   Adj, N, lgammaLookup, logLookup,
                                   DegreeCorrection,
                                  alpha, beta);

    //NOTE: When HotTemp == 0, ProbSwap is -NaN, and the swap
    //will never be accepted. This means that the random walk chain
    //never percolates through the other chains, so it basically just
    //wastes computing time. Should not have a large effect on results.
    ProbSwap = exp(((1/HotTemp) - (1/ColdTemp)) * (HotLik - ColdLik));

    if(SWAPS % 1000 == 0){
        fprintf(stderr, "percent swaps = %.4f\n", (double)ACCEPTEDSWAPS/(double)SWAPS);
    }

    if(gsl_ran_flat(r, 0 ,1) < ProbSwap){
        //swap the temperatures
        SWAPS++;
        ACCEPTEDSWAPS++;
        return 1;
    }
    else {
        SWAPS++;
        return 0;
    }
}
