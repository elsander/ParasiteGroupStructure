extern int PrintMatrixShort(FILE * F, gsl_matrix_short * A);
extern int PrintVectorFloat(FILE * F, gsl_vector * V);
extern int PrintVectorShort(FILE * F, gsl_vector_short * V);
extern int Partition_Initialize(gsl_vector_short * V,
                                int N,
                                int maxGroups,
                                gsl_rng * r);
extern int RGF(int N, gsl_vector_short * Chain, gsl_vector_short * RGFswap);
extern double Marginal(gsl_vector_short * V,
                       gsl_vector_short * VCopy,
                       gsl_matrix_short * Adj,
                       int N,
                       gsl_vector * lgammaLookup,
                       gsl_vector * logLookup);
extern double DegreeMarginal(gsl_vector_short * V,
                             gsl_vector_short * VCopy,
                             gsl_matrix_short * Adj,
                             int N,
                             double alpha,
                             double beta,
                             gsl_vector * lgammaLookup);
extern double Complexity(int X);
extern double NML(gsl_vector_short * V,
                  gsl_vector_short * VCopy,
                  gsl_matrix_short * Adj,
                  int N);
extern double Partition_Likelihood(gsl_vector_short * V,
                                   gsl_vector_short * VCopy,
                                   gsl_matrix_short * Adj,
                                   int N,
                                   gsl_vector * lgammaLookup,
                                   gsl_vector * logLookup,
                                   int DegreeCorrection,
                                   double alpha,
                                   double beta);
extern int TrySwap(int N,
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
                   double beta);
