// comtrait.c - community trait structure algs

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"
#include "stats.h"

//ComTraitMetric - calculate metrics of trait variation in each sample
void ComTraitMetric(sample S, traits T, int SwapMethod, int XVAR)
{

  //initialize variables
  int samp,trait, rec, run;
  int *STattach;
  int **rankLowSTMetricRnd;
  int **rankHighSTMetricRnd;
  float **STMeanObs;
  float ***STMeanRnd;
  float **MeanSTMetricRnd;
  float **STMetricObs;
  float ***STMetricRnd;
  float **StDevSTMetricRnd;
  float **SESSTMetricRnd;
  float *sampleTraits;
  float *randomMetrics;
  float mean, var, metric;

  rankLowSTMetricRnd = imatrix(0,S.nsamples-1,0,T.ntraits-1);
  rankHighSTMetricRnd = imatrix(0,S.nsamples-1,0,T.ntraits-1);

  STMeanObs = matrix(0,S.nsamples-1,0,T.ntraits-1);
  STMeanRnd = f3tensor(0,S.nsamples-1,0,T.ntraits-1,0,RUNS);
  MeanSTMetricRnd = matrix(0,S.nsamples-1,0,T.ntraits-1);

  STMetricObs = matrix(0,S.nsamples-1,0,T.ntraits-1);
  STMetricRnd = f3tensor(0,S.nsamples-1,0,T.ntraits-1,0,RUNS);
  StDevSTMetricRnd = matrix(0,S.nsamples-1,0,T.ntraits-1);

  SESSTMetricRnd = matrix(0,S.nsamples-1,0,T.ntraits-1);

  // attach sample to traits
  STattach = ivector(0, T.ntaxa-1);
  AttachSampleToTraits(S,T,STattach);

  //loop through samples and traits
  //calculate trait metric of each sample
  //store in an array
  for (trait = 0; trait < T.ntraits; trait++)
    {
      for (samp = 0; samp < S.nsamples; samp++)
        {
          //(re)dimension a vector to fill with trait values for this sample
          sampleTraits = vector(0,S.srec[samp]-1);

          //fill the vector
          for (rec = 0;rec < S.srec[samp];rec++)
            {
              sampleTraits[rec] = T.tr[STattach[ S.id[samp][rec] ]][trait];
            }

          //calculate the metric of the trait vector for this sample/trait combo
          meanvar(sampleTraits, S.srec[samp], &mean, &metric);
          STMeanObs[samp][trait] = mean;
          traitMetric(sampleTraits, S.srec[samp], &metric, XVAR);
          STMetricObs[samp][trait] = metric;
          free_vector(sampleTraits,0,S.srec[samp]-1);
        }
    }

  // if BURNIN > 0 and independent or trial swap null, randomize and discard (burn in)
  if ((BURNIN > 0) && (SwapMethod==3 || SwapMethod==4)) {
    //burn in using appropriate algorithm
    if (SwapMethod == 3) {
      IndependentSwap(S, BURNIN);
    } else {
      TrialSwap(S, BURNIN);
    }
  }

  //now repeat calculations for randomized matrices
  for (run = 0;run<RUNS;run++)
    {

      //Swap the sample or trait labels
      switch (SwapMethod)
        {
        case 0:
          //SwapMethod 0 = shuffle traits across species
          //(sample remains unshuffled)
          TraitsAttachShuffle(S,T,STattach);
          break;
        case 1:
          //SwapMethod 1 = samples become random draws from sample taxa list
          //(maintains sample species richnesses but not species frequencies)
          RandomizeSampleTaxaShuffle(S);
          break;
        case 2:
          //SwapMethod 2 = samples become random draws from traits
          //note this is redundant for traits but included for completeness
          RandomizeSampleTaxaShuffle(S);
          TraitsAttachShuffle(S,T,STattach);
          break;
        case 3:
          //SwapMethod 3 = independent checkerboard swap of sample matrix
          //(maintains species frequencies and sample species richnesses)
          IndependentSwap(S,SWAPS);
          break;
        case 4:
          //SwapMethod 4 = trial swap of sample matrix
          //(maintains species frequencies and sample species richnesses)
          TrialSwap(S, SWAPS);
          break;
        default:
          printf("Please use -m command line switch to specify a randomization method.\n");
          printf("See documentation for a list of possible null models.\n");
          exit(EXIT_FAILURE);
          break;
        }

      for (trait = 0; trait < T.ntraits; trait++)
        {
          for (samp = 0; samp < S.nsamples; samp++)
            {
              //(re)dimension a vector to fill with trait values for this sample
              sampleTraits = vector(0,S.srec[samp]-1);

              //fill the vector
              for (rec = 0;rec < S.srec[samp];rec++)
                {
                  sampleTraits[rec] = T.tr[STattach[ S.id[samp][rec] ]][trait];
                }

              //calculate the metric of the trait vector for this sample/trait combo
              meanvar(sampleTraits, S.srec[samp], &mean, &metric);
              STMeanRnd[samp][trait][run] = mean;
              traitMetric(sampleTraits, S.srec[samp], &metric, XVAR);
              STMetricRnd[samp][trait][run] = metric;
              free_vector(sampleTraits,0,S.srec[samp]-1);
            }
        }
    }

  //calculate summary statistics for randomized data
  for (trait = 0;trait < T.ntraits;trait++)
    {
      for (samp = 0;samp < S.nsamples;samp++)
        {
          randomMetrics = vector(0,RUNS-1);
          for (run = 0;run < RUNS;run++)
            {
              randomMetrics[run] = STMetricRnd[samp][trait][run];
              if (randomMetrics[run] > STMetricObs[samp][trait]) rankHighSTMetricRnd[samp][trait]++;
              if (randomMetrics[run] < STMetricObs[samp][trait]) rankLowSTMetricRnd[samp][trait]++;
            }
          meanvar(randomMetrics,RUNS,&mean,&var);
          MeanSTMetricRnd[samp][trait] = mean;
          StDevSTMetricRnd[samp][trait] = sqrt(var);
          SESSTMetricRnd[samp][trait] = (STMetricObs[samp][trait] - MeanSTMetricRnd[samp][trait]) \
            / StDevSTMetricRnd[samp][trait];
          free_vector(randomMetrics,0,RUNS-1);
        }
    }

  //summarize and output

  //print header
  printf("Phylocom output: randomization method %d, %d runs, trait metric %d (",SwapMethod,RUNS,XVAR);
  switch (XVAR)
    {
    case 1:
      printf("variance)\n");
      break;
    case 2:
      printf("MPD)\n");
      break;
    case 3:
      printf("MNTD)\n");
      break;
    case 4:
      printf("range)\n");
      break;
    default:
      printf(")\n");
      break;
    }

  //print results with separate entry for each trait/plot
  printf("Trait\tSample\tNTaxa\tMean\tMetric\tMeanRndMetric\tSDRndMetric\tSESMetric\trankLow\trankHigh\truns\n");
  for (trait = 0;trait < T.ntraits;trait++)
    {
      for (samp = 0;samp < S.nsamples;samp++)
        {
          printf("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\n",T.trname[trait],S.pname[samp],S.srec[samp],STMeanObs[samp][trait],STMetricObs[samp][trait],MeanSTMetricRnd[samp][trait],StDevSTMetricRnd[samp][trait],SESSTMetricRnd[samp][trait],rankLowSTMetricRnd[samp][trait],rankHighSTMetricRnd[samp][trait],RUNS);
        }
    }

}

// --------- Phylogeny Attach Shuffle ---------------------------
void TraitsAttachShuffle(sample S, traits T, int *attach)
{

  //This randomization works by shuffling the attach vector connecting
  //sample taxa to their correponding trait values
  //This effectively randomizes trait values across species

  int taxon;
  int randomTraitsTaxon;
  int oldTaxon;

  //assign each taxon in sample a new random corresponding trait taxon#
  for (taxon = 0;taxon<S.ntaxa;taxon++)
    {
      randomTraitsTaxon = ((int)(((double)T.ntaxa*random())/(double)(RAND_MAX)));
      oldTaxon = attach[taxon];
      attach[taxon] = randomTraitsTaxon;
      attach[randomTraitsTaxon] = oldTaxon;
    }
}

void traitMetric(float X[], unsigned int n, float *metric, int XVAR)
{
  //calculate various metrics of trait dispersion in samples
  //set metric with the -x command line switch (fix later)
  //default XVAR = 0 = variance of traits in each sample
  unsigned long i,j;
  float minTrait, maxTrait;
  long double sumx = 0.0; // sum of X - long, just in case
  long double sumxsq = 0.0; // sum 0f X^2
  long double totalDist = 0.0; //distance counter
  double minDist = 0.0; //minimum distance
  unsigned long numPairs = 0; //pairs counter

  switch (XVAR)
    {
    case 1:
      //Trait Metric 0 = calculate variance of traits in each sample
      for (i = 0; i < n; i++)
        {
          sumx   += X[i];
          sumxsq += X[i] * X[i];
        }
      *metric = (float) (( sumxsq - ( (sumx * sumx) / (double) n ) ) / (double) (n-1));
      break;
    case 2:
      //Trait Metric 1 = calculate mean pairwise distance among traits in each sample
      for (i = 0; i < n - 1; i++)
        {
          for (j = 1; j < n; j++)
            {
              totalDist += fabs(X[i] - X[j]);
              numPairs++;
            }
        }
      *metric = (float) (totalDist / (double) numPairs);
      break;
    case 3:
      //Trait Metric 2 = calculate mean nearest neighb. dist among traits in samples
      for (i = 0; i < n; i++)
        minDist = 99999999.9; //needs to be higher than minimum expected distance
      for (j = 0; j < n; j++)
        {
          if (i != j)
            {
              //find lowest pairwise distance between taxa
              if (fabs(X[i] - X[j]) < minDist) minDist = fabs(X[i] - X[j]);
            }
          totalDist += minDist;
        }
      *metric = (float) ( totalDist / (double) n );
      break;
    case 4:
      //Trait Metric 3 = calculate trait range in sample
      minTrait = X[0];
      maxTrait = X[0];
      for (i = 0; i < n; i++)
        {
          if (X[i] < minTrait) minTrait = X[i];
          if (X[i] > maxTrait) maxTrait = X[i];
        }
      *metric = (float) (maxTrait - minTrait);
      break;
    default:
      printf("Please use -x command line switch to specify a trait metric.\n");
      printf("See documentation for a list of possible metrics.\n");
      exit(EXIT_FAILURE);
      break;
    }
}
