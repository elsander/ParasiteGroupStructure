#include "Common.h"
#include "partition.h"

//hardcoded variables
//Note that HOTTEMP should be greater than 0
#define HOTTEMP .0000001 //hottest temperature for MC3
#define COLDTEMP 1 //coldest temperature for MC3
#define SWAPSTEPS 20 //how often to try to swap temperatures

//Gibbs sampler step
double Gibbs (int N,
              gsl_vector_short * Chain,
              gsl_vector_short * ChainCopy,
              gsl_vector_short * ChainCopy2,
              gsl_matrix_short * Adj,
              double Temp,
              gsl_vector_short * RGFswap,
              gsl_rng * r,
              gsl_vector * lgammaLookup,
              gsl_vector * logLookup,
              int maxGroups,
              int DegreeCorrection,
              double alpha,
              double beta
    ){

    //add one to nGroups to allow for a new group to be formed
    int nGroups = gsl_vector_short_max(Chain);
	//add one unless it already has maxGroups
    if(nGroups < maxGroups) {
        nGroups++;
    }

    //initialize vector of swap likelihoods
    //it will need to be length maxGroups at most
    double swapLikelihoods[maxGroups];

//choose a random vector element to mutate
    int mutInd = gsl_rng_uniform_int(r, N);
    int i;
    //get likelihoods for all possible swaps for that element
    //this includes no change
    for(i=0; i<=nGroups; i++){
        gsl_vector_short_memcpy(ChainCopy, Chain);
        gsl_vector_short_set(ChainCopy, mutInd, i);
        RGF(N, ChainCopy, RGFswap);
        swapLikelihoods[i] = Temp * Partition_Likelihood(ChainCopy,
                                                         ChainCopy2,
                                                         Adj,
                                                         N,
                                                         lgammaLookup,
                                                         logLookup,
                                                         DegreeCorrection,
                                                         alpha,
                                                         beta);
    }

    //find the minimum in the array
    double minLik = 1000000;
    for(i=0; i<nGroups; i++){
        if(swapLikelihoods[i] < minLik){
            minLik = swapLikelihoods[i];
        }
    }

    double CutOff = .0000000000000001; //precision cutoff = 10e-16
    CutOff = log(CutOff) - log(nGroups);

    //subtract largest value from all, to avoid precision issues
    for(i=0; i<nGroups; i++){
        swapLikelihoods[i] = swapLikelihoods[i] - minLik;
        if(swapLikelihoods[i] < CutOff){
            swapLikelihoods[i] = sqrt(-1); //should result in NaN
        }
        swapLikelihoods[i] = exp(swapLikelihoods[i]);
    }

    //turn these likelihoods into cumulative probabilities for
    //sampling
    double runningTotal = 0.0;
    for(i=0; i<nGroups; i++){
        if(!isnan(swapLikelihoods[i])){
            runningTotal += swapLikelihoods[i];
            swapLikelihoods[i] = runningTotal;
        }
    }

    //choose a mutation based on these weighted probabilities
    double mutVal = gsl_ran_flat(r, 0.0, runningTotal);
    //figure out which mutation this corresponds to
    int chosenIndex = 0;
    for(i=0; i<nGroups; i++){
        if(mutVal <= swapLikelihoods[i]){
            chosenIndex = i;
            break;
        }
    }

    //now update the Chain
    gsl_vector_short_set(Chain, mutInd, chosenIndex);
    RGF(N, Chain, RGFswap);

    return Partition_Likelihood(Chain, ChainCopy, Adj,
                                N, lgammaLookup, logLookup,
                                DegreeCorrection,
                                alpha, beta);
}

//MCMCMC algorithm
double MC3 (int N,
            gsl_matrix_short * Adj,
            int Steps,
            int nChains,
            gsl_vector_short * BestSolution,
            gsl_rng * r,
            gsl_vector * lgammaLookup,
            gsl_vector * logLookup,
            int maxGroups,
            int DegreeCorrection,
            double alpha,
            double beta
    ){

    // create the chains
    gsl_vector_short * Chains[nChains];
    //create copies for use by Gibbs and marginal functions
    gsl_vector_short * ChainCopy = gsl_vector_short_calloc(N);
    gsl_vector_short * ChainCopy2 = gsl_vector_short_calloc(N);
    // create the fitness vector
    gsl_vector * Marginals = gsl_vector_calloc(nChains);
    //initialize swapping vector for RGF
    gsl_vector_short * RGFswap = gsl_vector_short_calloc(N+1);

    int i,j,k;
    double BestMarginal;
    BestMarginal = -1000000000.0;

    //initialize chains
    for(i=0; i<nChains; i++){
        // allocate memory
        Chains[i] = gsl_vector_short_calloc(N);
        // initialize the population
        Partition_Initialize(Chains[i], N, maxGroups, r);
        RGF(N, Chains[i], RGFswap);
    }

    //generate temperatures assuming uniform spacing
    gsl_vector * Temps = gsl_vector_calloc(nChains);
    //step size for incrementing temperatures
    double StepSize;
    StepSize = (COLDTEMP - HOTTEMP)/((double)nChains - 1);
    gsl_vector_set(Temps, 0, HOTTEMP);
    for(i=1; i<(nChains-1); i++){
        gsl_vector_set(Temps, i, gsl_vector_get(Temps, i-1)+StepSize);
    }
    gsl_vector_set(Temps, nChains-1, COLDTEMP);

    //for convenience, we want a copy of the Temps vector that doesn't
    //get swapped around
    gsl_vector * TempsCopy = gsl_vector_calloc(nChains);
    gsl_vector_memcpy(TempsCopy, Temps);

    int chInd1, chInd2;
    double dtmp;
    int itmp;
    int swapFlag;

    for(i=0; i<Steps; i++){
        //print the best likelihood we've found so far every so often
        if(i % 1000 == 0){
            fprintf(stderr, "Step %d Best solution %1.4f\n", i, BestMarginal);
        }

        //if enough steps have passed, swap temperatures
        if(i % SWAPSTEPS == 0){
            //try to swap using "bucket brigade"
            for(j=0; j<(nChains-1); j++){
                //find which chains have adjacent temperatures
                chInd1 = -1;
                chInd2 = -1;
                for(k=0; k<nChains; k++){
                    if(gsl_vector_get(TempsCopy, j) == gsl_vector_get(Temps, k)){
                        chInd1 = k;
                    }
                    if(gsl_vector_get(TempsCopy, j+1) == gsl_vector_get(Temps, k)){
                        chInd2 = k;
                    }
                    if(chInd1 >= 0 && chInd2 >=0){
                        break;
                    }
                }
                //try to swap them
                swapFlag = TrySwap(N, Adj,
                                   Chains[chInd1], Chains[chInd2],
                                   ChainCopy,
                                   gsl_vector_get(Temps, chInd1),
                                   gsl_vector_get(Temps, chInd2),
                                   r, lgammaLookup, logLookup,
                                   DegreeCorrection,
                                   alpha, beta);
                if(swapFlag == 1){
                    dtmp = gsl_vector_get(Temps, chInd1);
                    gsl_vector_set(Temps, chInd1, gsl_vector_get(Temps, chInd2));
                    gsl_vector_set(Temps, chInd2, dtmp);
                }
            }
        }

        //take a step
        for(j=0; j<nChains; j++){
            dtmp = Gibbs(N, Chains[j], ChainCopy,
                         ChainCopy2, Adj,
                         gsl_vector_get(Temps, j), RGFswap,
                         r, lgammaLookup, logLookup, maxGroups,
                         DegreeCorrection, alpha, beta);
            gsl_vector_set(Marginals, j, dtmp);
        }

        //update the best solution, if appropriate
        if(gsl_vector_max(Marginals) > BestMarginal){
            itmp = gsl_vector_max_index(Marginals);
            BestMarginal = gsl_vector_get(Marginals, itmp);
            gsl_vector_short_memcpy(BestSolution, Chains[itmp]);
            fprintf(stderr, "Steps %d Best solution %.4f\n", i, BestMarginal);
        }
    }

    //free memory
    gsl_vector_short_free(RGFswap);
    gsl_vector_free(Temps);
    gsl_vector_free(TempsCopy);
    gsl_vector_free(Marginals);
    gsl_vector_short_free(ChainCopy);
    gsl_vector_short_free(ChainCopy2);
    for(i=0; i<nChains; i++){
        gsl_vector_short_free(Chains[i]);
    }

    return BestMarginal;
}
