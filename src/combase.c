// combase.c - main community phylogenetic structure algorithms

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

// ---------- RELATEDNESS -----------------------------------------------------
float Relatedness(struct phylo P, int phyTaxonA, int phyTaxonB)
// A,B are phylo taxalist numbers (internal and terminal)
// fixed cw-11jun05
{
  float cumulativeDistance = 0.0;
  int nodeA, nodeB;

  nodeA = P.t2n[phyTaxonA];
  nodeB = P.t2n[phyTaxonB];

  while (nodeA != nodeB) {
    if (P.depth[nodeA] >= P.depth[nodeB]) {
      cumulativeDistance += P.bl[nodeA];
      nodeA = P.up[nodeA];
    } else {
      cumulativeDistance += P.bl[nodeB];
      nodeB = P.up[nodeB];
    }
  }

  return cumulativeDistance;
}

// ---------- TaxaDist ------------------------------------------------------
// Simplifies lookup of phylogenetic distances given sample taxa codes
float TaxaDist(struct phylo P, int sampTaxonA, int sampTaxonB, int *attach) {
  return P.dist[P.t2n[attach[sampTaxonA]]][P.t2n[attach[sampTaxonB]]];
}

// ---------- SIMPLEDIST ------------------------------------------------------
void SimpleDist(phylo Intree) {
  // t2n fixed cw

  int i, j;

  // allocate a matrix of memory to the pre-existing pointer phylo.dist
  Intree.dist = matrix(0, Intree.nnodes-1, 0, Intree.nnodes-1);
  DistMatrix(Intree);

  // headers
  printf(".");
  for (i= 0; i < Intree.nnodes; i++) {
    if (Intree.noat[i] == 0)
      printf("\t%s", Intree.taxon[i]);
  }
  printf("\n");

  // rows
  for (i= 0; i < Intree.nnodes; i++) {
    if (Intree.noat[i] == 0) {
      //print row name
      printf("%s", Intree.taxon[i]);

      // print distances in columns
      for (j= 0; j < Intree.nnodes; j++) {
        if (Intree.noat[j] == 0)
          printf("\t%f", Intree.dist[i][j]);
      }
      printf("\n");
    }
  }
}

// ---------- PHYLOVARCOVAR ------------------------------------------------------
void PhyloVarCovar(phylo Intree) {
  // Calculates phylogenetic variance-covariance matrix
  // Results output as square matrix
  // Diagonals - variance - branch length from root node to taxon
  // Off-diagonals - covariance - branch length from root to tips shared by 2 taxa
  // This function is SK's modification of CW's SimpleDist() function

  int i, j, MRCA;

  // headers
  printf(".");
  for (i= 0; i < Intree.nnodes; i++) {
    if (Intree.noat[i] == 0)
      printf("\t%s", Intree.taxon[i]);
  }
  printf("\n");

  // rows
  for (i= 0; i < Intree.nnodes; i++) {
    if (Intree.noat[i] == 0) {
      //print row name
      printf("%s", Intree.taxon[i]);

      // print distances in columns
      for (j= 0; j < Intree.nnodes; j++) {
        if (Intree.noat[j] == 0) {
          if (i == j) {
            //diagonals - "variance" = distance from taxon to root node
            printf("\t%f", DistToRootNode(Intree, i));
          } else {
            //off-diagonals - "covariance" = dist from taxa's MRCA to root node
            //find MRCA for two nodes
            MRCA = FindMRCA(Intree, i, j);
            printf("\t%f", DistToRootNode(Intree, MRCA));
          }
        }
      }
      printf("\n");
    }
  }
}

// ---------- ATTACHSAMPLETOPHYLO ----------------------------------------
void AttachSampleToPhylo(struct sample S, struct phylo P, int *SPattach) {
  //creates lookup table
  //SPattach[sample.taxa#] returns phylo taxalist # of that sample taxon
  int i, j;
  int found = 0;

  for (i = 0; i < S.ntaxa; i++) {
    found = 0;
    for (j = 0; j < P.termtaxa; j++) {
      if (strcmp(S.taxa[i], P.taxalist[j]) == 0) {
        SPattach[i] = j;
        found = 1;
      }
    }
    if (found == 0) {
      printf(
             "Taxon %s in sample not a terminal taxon in phylogeny. Exiting.\n",
             S.taxa[i]);
      exit(EXIT_FAILURE);
    }
    if (Debug)
      printf("Sid %d\tPid %d\tSn %s\tPn %s\n", i, SPattach[i], S.taxa[i], P.taxalist[SPattach[i]]);
  }

}

// ---------- ATTACH ---------------------------------------------------------
void AttachSampleToTraits(struct sample S, struct traits T, int *STattach) {
  //creates lookup table - STattach[sample.taxa#] = traits.taxon #
  int i, j;
  int found = 0;

  for (i = 0; i < S.ntaxa; i++) {
    found = 0;
    for (j = 0; j < T.ntaxa; j++) {
      if (strcmp(S.taxa[i], T.taxon[j]) == 0) {
        STattach[i] = j;
        found = 1;
      }
    }
    if (found == 0) {
      printf("Taxon %s in sample not found in traits. Exiting.\n", S.taxa[i]);
      exit(EXIT_FAILURE);
    }
  }

}

// ---------- ATTACHTRAITSTOPHYLO ----------------------------------------
void AttachTraitsToPhylo(struct traits T, struct phylo P, int *TPattach) {
  //creates lookup table - TPattach[traits.taxon#] = phylo.taxalist #
  int i, j;
  int found = 0;

  for (i = 0; i < T.ntaxa; i++) {
    found = 0;
    for (j = 0; j < P.termtaxa; j++) {
      if (strcmp(T.taxon[i], P.taxalist[j]) == 0) {
        TPattach[i] = j;
        found = 1;
      }
    }
    if (found == 0) {
      printf("Taxon %s in traits not found in phylo. Exiting.\n", T.taxon[i]);
      exit(EXIT_FAILURE);
    }
  }
}

// ------------ ATTACHPHYLOTOTRAITS ----------------------
void AttachPhyloToTraits(struct phylo P, struct traits T, int *PTattach) {
  //lookup table - PTattach[phylo.node#] returns Trait.taxon# of [node] in phylo
  int i, j, matched;

  for (i = 0; i < T.ntaxa; i++) {
    matched = 0;
    // cw fixed ntaxa
    for (j = 0; j < P.ntaxa; j++) {
      if (strcmp(P.taxalist[j], T.taxon[i]) == 0) {
        PTattach[P.t2n[j]] = i;
        matched = 1;
      }
    }

    if (matched == 0) {
      printf(
             "  Taxa string `%s' in traits file not found in phylo file.\n",
             T.taxon[i]);
      printf("  Exiting.\n");
      exit(EXIT_FAILURE);
    }
  }
}

//---------- DistMatrix -------------
void DistMatrix(struct phylo P) {
  // fills .dist matrix of all node-node distances for phylo P
  // i and j are nodes in the phylogeny
  //TODO check to make sure works with non-ultrametric trees
  int i, j;
  float cumulativeDist;
  int xnode, ynode;

  for (i = 0; i < P.nnodes; i++) {
    for (j = 0; j < P.nnodes; j++) {
      cumulativeDist = 0.0;
      xnode = i;
      ynode = j;

      while (xnode != ynode) {
        if (P.depth[xnode] >= P.depth[ynode]) {
          cumulativeDist += P.bl[xnode];
          xnode = P.up[xnode];
        } else {
          cumulativeDist += P.bl[ynode];
          ynode = P.up[ynode];
        }
      }

      P.dist[i][j] = cumulativeDist;
      //printf("%d,%d\t%f\n", i, j, cumulativeDist);
    }
  }
}

//---------- DistMatrixNN -------------
void DistMatrixNN(struct phylo P) {
  // fills .dist matrix of all node-node distances for phylo P
  // i and j are nodes in the phylogeny
  // SWK's modification of CW's DistMatrix()
  // Modified to work with non-terminal nodes and non-ultrametric trees

  int i, j;
  float cumulativeDist;
  int xNode, yNode;
  int xDepth, yDepth;

  for (i = 0; i < P.nnodes; i++) {
    for (j = 0; j < P.nnodes; j++) {
      cumulativeDist = 0.0;
      xNode = i;
      yNode = j;

      while (xNode != yNode) {
        xDepth = P.depth[xNode];
        yDepth = P.depth[yNode];
        if (xDepth >= yDepth) {
          cumulativeDist += P.bl[xNode];
          xNode = P.up[xNode];
        }
        if (yDepth >= xDepth) {
          cumulativeDist += P.bl[yNode];
          yNode = P.up[yNode];
        }
      }

      P.dist[i][j] = cumulativeDist;
      if (Debug)
        printf("nodes %d %s,%d %s:\tdistance %f\n", i, P.taxon[i], j,
               P.taxon[j], cumulativeDist);
    }
  }
}

// ---------- DIST TO ROOT NODE -----------
// Calculates distance (depth in branch length units) from node x to root node
float DistToRootNode(phylo P, int node) {
  float cumulativeDist = 0.0;

  while (node != 0) {
    cumulativeDist += P.bl[node];
    node = P.up[node];
  }

  return cumulativeDist;
}

// ---- FIND MOST RECENT COMMON ANCESTOR ----
// Finds node that is the MRCA of two other nodes in the phylo
int FindMRCA(phylo P, int xNode, int yNode) {
  int xDepth, yDepth;

  //Move up through phylo until at deepest common node - this is MRCA
  while (xNode != yNode) {
    xDepth = P.depth[xNode];
    yDepth = P.depth[yNode];
    if (xDepth >= yDepth)
      xNode = P.up[xNode];
    if (yDepth >= xDepth)
      yNode = P.up[yNode];
  }

  return xNode;
}

// ----------- NEW CLUST --------------------------------------------------
void ComStruct(phylo InP, sample InS, int SwapMethod, int UseAbund) {
  //This function written by Steven Kembel
  //Calculates community phylogenetic structure
  //Tests significance against null communities generated
  //using various swapping algorithms

  //Now works using one of several SwapMethods.

  int plot;
  int runNum;
  double NR_plot;
  double NT_plot;
  int *pLowCountMeanX;
  int *pHighCountMeanX;
  int *pLowCountMeanXMin;
  int *pHighCountMeanXMin;
  double *meanX;
  double *randomized_meanX; //mean x in a plot across runs
  double *randomized_sdMeanX;
  double *meanXMin;
  double *randomized_meanXMin;
  double *randomized_sdMeanXMin;
  double *plot_NRI;
  double *plot_NTI;
  double **randomized_plotMeanX;
  double **randomized_plotMeanXMin;
  int *attach;

  //allocate vectors and matrices
  pLowCountMeanX = ivector(0, InS.nsamples-1);
  pHighCountMeanX = ivector(0, InS.nsamples-1);
  pLowCountMeanXMin = ivector(0, InS.nsamples-1);
  pHighCountMeanXMin = ivector(0, InS.nsamples-1);
  meanX = dvector(0, InS.nsamples-1);
  randomized_meanX = dvector(0, InS.nsamples-1); //mean x in a plot across runs
  randomized_sdMeanX = dvector(0, InS.nsamples-1);
  meanXMin = dvector(0, InS.nsamples-1);
  randomized_meanXMin = dvector(0, InS.nsamples-1);
  randomized_sdMeanXMin = dvector(0, InS.nsamples-1);
  plot_NRI = dvector(0, InS.nsamples-1);
  plot_NTI = dvector(0, InS.nsamples-1);
  randomized_plotMeanX = dmatrix(0, InS.nsamples-1, 0, RUNS);
  randomized_plotMeanXMin = dmatrix(0, InS.nsamples-1, 0, RUNS);

  //Attach sample to phylogeny
  attach = ivector(0, InS.ntaxa-1);
  AttachSampleToPhylo(InS, InP, attach);

  //Create and fill matrix of distances among all nodes in the phylogeny
  InP.dist = matrix(0, InP.nnodes-1, 0, InP.nnodes-1);
  DistMatrix(InP);

  //initialize counters
  for (plot = 0; plot < InS.nsamples; plot++) {
    randomized_meanX[plot] = 0.0;
    randomized_meanXMin[plot] = 0.0;
    randomized_sdMeanX[plot] = 0.0;
    randomized_sdMeanXMin[plot] = 0.0;
    pLowCountMeanX[plot]= 0;
    pHighCountMeanX[plot]= 0;
    pLowCountMeanXMin[plot] = 0;
    pHighCountMeanXMin[plot] = 0;
  }

  //output header line
  if (SwapMethod == 3 || SwapMethod == 4) {
    printf("Phylocom output: randomization method %d, %d swaps/trials,  %ld burnin swaps/trials, %d runs\n",
           SwapMethod, SWAPS, BURNIN, RUNS);
  } else {
    printf("Phylocom output: randomization method %d, %d runs\n",
           SwapMethod, RUNS);
  }

  //  printf("plot       ntaxa       MPD       NRI    pMPD      MNTD       NTI   pMNTD\n");

  // Calculate results for real data:
  for (plot = 0; plot < InS.nsamples; plot++) {
    //calc means for each plot
    meanX[plot] = MeanDistance(InP, InS, attach, plot, UseAbund);
    meanXMin[plot] = MeanMinimumDistance(InP, InS, attach, plot, UseAbund);
  }

  // if Verbose, output header line for verbose output
  if (Verbose)
    printf("Run\tSample\tMPD\tMNTD\n");

  // if BURNIN > 0 and independent or trial swap null, randomize and discard (burn in)
  if ((BURNIN > 0) && (SwapMethod==3 || SwapMethod==4)) {
    //burn in using appropriate algorithm
    if (SwapMethod == 3) {
      IndependentSwap(InS, BURNIN);
    } else {
      TrialSwap(InS, BURNIN);
    }
  }

  // Now get results for a number of randomized null communities
  for (runNum = 0; runNum < RUNS; runNum++) {
    //Swap the matrix or phylogeny
    switch (SwapMethod) {
    case 0:
      //SwapMethod 0 = shuffle taxa labels on phylogeny
      //(sample remains unshuffled)
      PhylogenyAttachShuffle(InP, InS, attach);
      break;
    case 1:
      //SwapMethod 1 = samples become random draws from sample taxa list
      //(maintains sample species richnesses but not species frequencies)
      //(Species present in the phylogeny but not the sample will not
      //be represented in the shuffled samples)
      RandomizeSampleTaxaShuffle(InS);
      break;
    case 2:
      //SwapMethod 2 = samples become random draws from phylogeny
      //(maintains sample species richnesses but not species frequencies)
      //(also shuffles link between phylo and sample taxa to give
      //species present in phylo but not sample chance to be represented)
      //Note limitation of this way of shuffling is that if phylogeny
      //contains more taxa than sample, each shuffle will only
      //fill the sample with sample # taxa, not phylo # taxa
      RandomizeSampleTaxaShuffle(InS);
      PhylogenyAttachShuffle(InP, InS, attach);
      break;
    case 3:
      //SwapMethod 3 = independent checkerboard swap of sample matrix
      //(maintains species frequencies and sample species richnesses)
      IndependentSwap(InS, SWAPS);
      break;
    case 4:
      //SwapMethod 4 = trial swap of sample matrix
      //(maintains species frequencies and sample species richnesses)
      TrialSwap(InS, SWAPS);
      break;
    default:
      printf("Please use -m command line switch to specify a randomization method.\n");
      printf("See documentation for a list of possible null models.\n");
      exit(EXIT_FAILURE);
      break;
    }

    // Get results for randomized data
    for (plot = 0; plot < InS.nsamples; plot++) {
      //calculate NRI- and NTI-type metrics for the plot

      NR_plot = MeanDistance(InP, InS, attach, plot, UseAbund);
      randomized_plotMeanX[plot][runNum] = NR_plot;
      randomized_meanX[plot] += NR_plot;

      if (NR_plot >= meanX[plot]) {
        pLowCountMeanX[plot]++;
      }
      if (NR_plot <= meanX[plot]) {
        pHighCountMeanX[plot]++;
      }

      NT_plot = MeanMinimumDistance(InP, InS, attach, plot, UseAbund);
      randomized_plotMeanXMin[plot][runNum] = NT_plot;
      randomized_meanXMin[plot] += NT_plot;
      if (NT_plot >= meanXMin[plot]) {
        pLowCountMeanXMin[plot]++;
      }
      if (NT_plot <= meanXMin[plot]) {
        pHighCountMeanXMin[plot]++;
      }

      //if Verbose output desired, write raw random MPD/MNTD
      if (Verbose)
        printf("%d\t%s\t%f\t%f\n", runNum, InS.pname[plot], NR_plot,
               NT_plot);

    }
  }

  // Calculate mean and standard deviation of randomized meanX/meanXMin for each plot
  for (plot = 0; plot < InS.nsamples; plot++) {
    randomized_meanX[plot] /= (double)(RUNS);
    randomized_meanXMin[plot] /= (double)(RUNS);
    //Tally sums of squares - divide by N-1 during tally to avoid overflows
    for (runNum = 0; runNum < RUNS; runNum++) {
      randomized_sdMeanX[plot] += pow(randomized_meanX[plot]
                                      -randomized_plotMeanX[plot][runNum], 2)/ (double)(RUNS - 1);
      randomized_sdMeanXMin[plot] += pow(randomized_meanXMin[plot]
                                         -randomized_plotMeanXMin[plot][runNum], 2)/ (double)(RUNS
                                                                                              - 1);
    }
    randomized_sdMeanX[plot] = sqrt(randomized_sdMeanX[plot]);
    randomized_sdMeanXMin[plot] = sqrt(randomized_sdMeanXMin[plot]);
  }

  // Output header line
  printf("plot\tntaxa\tMPD\tMPD.rnd\tMPD.sd\tNRI\tMPD.rankLow\tMPD.rankHi\tMNTD\tMNTD.rnd\tMNTD.sd\tNTI\tMNTD.rankLo\tMNTD.rankHi\truns\n");

  // Calculate low/high p-value for meanX and meanXMin (for plots and total)
  // (see Legendre & Legendre 1998)
  for (plot = 0; plot < InS.nsamples; plot++) {

    //Calculate plot NRI/NTI
    plot_NRI[plot] = (double)(-1.0) * ( (meanX[plot]
                                         - randomized_meanX[plot]) / randomized_sdMeanX[plot]);
    plot_NTI[plot] = (double)(-1.0) * ( (meanXMin[plot]
                                         - randomized_meanXMin[plot]) / randomized_sdMeanXMin[plot]);

    //Output results tab-delimited
    printf(
           "%s\t%i\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%i\t%i\t%3.4f\t%3.4f\t%3.4f\t%3.4f\t%i\t%i\t%i\n",
           InS.pname[plot], InS.srec[plot], meanX[plot],
           randomized_meanX[plot], randomized_sdMeanX[plot],
           plot_NRI[plot], pLowCountMeanX[plot], pHighCountMeanX[plot],
           meanXMin[plot], randomized_meanXMin[plot],
           randomized_sdMeanXMin[plot], plot_NTI[plot],
           pLowCountMeanXMin[plot], pHighCountMeanXMin[plot], RUNS);
  }
}

// ---------- INDEPENDENT SWAP ----------------------------------------------------

void IndependentSwap(sample S, int NumberOfSwaps) {
  int swapped;
  int swapNum;
  int randomTaxon1, randomTaxon2, randomPlot1, randomPlot2;
  int idHolder;
  int abundHolder;
  int SppOnePresentInPlotOne, SppOnePresentInPlotTwo;
  int SppTwoPresentInPlotOne, SppTwoPresentInPlotTwo;
  int individualOne, individualTwo;
  int sppOneIndiv = 0, sppOnePlot = 0;
  int sppTwoIndiv = 0, sppTwoPlot = 0;

  //this procedure will shuffle the raw data using independent swap algorithm
  //holding plot species richness and species relative abundance constant
  //(the "independent swap" described in Gotelli and Entsminger)
  //(a.k.a. SIM9 in Gotelli (2000) Ecology)

  for (swapNum = 0; swapNum < NumberOfSwaps; swapNum++) {
    swapped = FALSE;

    //TODO add option to do the "Trial Swap" method (Miklos and Podani 2004)?

    //numTries = 0;
    while (swapped == FALSE) {
      //try to find and swap a checkerboard pattern in the matrix

      //Pick 2 random plots and 2 random taxa
      randomTaxon1 = random() % S.ntaxa;
      randomTaxon2 = random() % S.ntaxa;
      while (randomTaxon1 == randomTaxon2) {
        randomTaxon2 = ((int) (((double) S.ntaxa * random())
                               / (double)(RAND_MAX)) );
      }
      randomPlot1 = random() % S.nsamples;
      randomPlot2 = random() % S.nsamples;
      while (randomPlot1 == randomPlot2) {
        randomPlot2 = random() % S.nsamples;
      }

      //check for checkerboard pattern of type
      //(0..1)  (1..0)
      //or
      //(1..0)  (0..1)
      //if it is a checkerboad, swap the values in Id and exit
      SppOnePresentInPlotOne = FALSE;
      SppOnePresentInPlotTwo = FALSE;
      SppTwoPresentInPlotOne = FALSE;
      SppTwoPresentInPlotTwo = FALSE;

      // count # of presence/absences in the 4 cells
      for (individualOne = 0; individualOne < S.srec[randomPlot1]; individualOne++) {
        for (individualTwo = 0; individualTwo < S.srec[randomPlot2]; individualTwo++) {
          if (S.id[randomPlot1][individualOne] == randomTaxon1) {
            SppOnePresentInPlotOne = TRUE;
            sppOneIndiv = individualOne;
            sppOnePlot = randomPlot1;
          }
          if (S.id[randomPlot2][individualTwo] == randomTaxon1) {
            SppOnePresentInPlotTwo = TRUE;
            sppOneIndiv = individualTwo;
            sppOnePlot = randomPlot2;
          }
          if (S.id[randomPlot1][individualOne] == randomTaxon2) {
            SppTwoPresentInPlotOne = TRUE;
            sppTwoIndiv = individualOne;
            sppTwoPlot = randomPlot1;
          }
          if (S.id[randomPlot2][individualTwo] == randomTaxon2) {
            SppTwoPresentInPlotTwo = TRUE;
            sppTwoIndiv = individualTwo;
            sppTwoPlot = randomPlot2;
          }
        }
      }

      //if checkboard, swap species between plots
      if (SppOnePresentInPlotOne == SppTwoPresentInPlotTwo) {
        if (SppOnePresentInPlotTwo == SppTwoPresentInPlotOne) {
          if (SppOnePresentInPlotOne != SppTwoPresentInPlotOne) {
            //if it is checkerboard, swap Id and Abund between plots
            //after storing current values
            idHolder = S.id[sppOnePlot][sppOneIndiv];
            S.id[sppOnePlot][sppOneIndiv]
              = S.id[sppTwoPlot][sppTwoIndiv];
            S.id[sppTwoPlot][sppTwoIndiv] = idHolder;
            abundHolder = S.abund[sppOnePlot][sppOneIndiv];
            S.abund[sppOnePlot][sppOneIndiv]
              = S.abund[sppTwoPlot][sppTwoIndiv];
            S.abund[sppTwoPlot][sppTwoIndiv] = abundHolder;
            swapped = TRUE;
          }
        }
      }

    }
  }
}


// ---------- TRIAL SWAP ----------------------------------------------------

void TrialSwap(sample S, int NumberOfTrials) {

  int trialNum;
  int randomTaxon1, randomTaxon2, randomPlot1, randomPlot2;
  int idHolder;
  int abundHolder;
  int SppOnePresentInPlotOne, SppOnePresentInPlotTwo;
  int SppTwoPresentInPlotOne, SppTwoPresentInPlotTwo;
  int individualOne, individualTwo;
  int sppOneIndiv = 0, sppOnePlot = 0;
  int sppTwoIndiv = 0, sppTwoPlot = 0;

  //this procedure will shuffle the raw data using trial swap algorithm
  //holding plot species richness and species relative abundance constant
  //(described by Miklos and Podani 2004)

  for (trialNum = 0; trialNum < NumberOfTrials; trialNum++) {

    //each trial: try to find and swap a checkerboard pattern in the matrix

    //Pick 2 random plots and 2 random taxa
    randomTaxon1 = random() % S.ntaxa;
    randomTaxon2 = random() % S.ntaxa;
    while (randomTaxon1 == randomTaxon2) {
      randomTaxon2 = ((int) (((double) S.ntaxa * random())
                             / (double)(RAND_MAX)) );
    }
    randomPlot1 = random() % S.nsamples;
    randomPlot2 = random() % S.nsamples;
    while (randomPlot1 == randomPlot2) {
      randomPlot2 = random() % S.nsamples;
    }

    //check for checkerboard pattern of type
    //(0..1)    (1..0)
    //or
    //(1..0)    (0..1)
    //if it is a checkerboad, swap the values in Id and exit
    SppOnePresentInPlotOne = FALSE;
    SppOnePresentInPlotTwo = FALSE;
    SppTwoPresentInPlotOne = FALSE;
    SppTwoPresentInPlotTwo = FALSE;

    // count # of presence/absences in the 4 cells
    for (individualOne = 0; individualOne < S.srec[randomPlot1]; individualOne++) {
      for (individualTwo = 0; individualTwo < S.srec[randomPlot2]; individualTwo++) {
        if (S.id[randomPlot1][individualOne] == randomTaxon1) {
          SppOnePresentInPlotOne = TRUE;
          sppOneIndiv = individualOne;
          sppOnePlot = randomPlot1;
        }
        if (S.id[randomPlot2][individualTwo] == randomTaxon1) {
          SppOnePresentInPlotTwo = TRUE;
          sppOneIndiv = individualTwo;
          sppOnePlot = randomPlot2;
        }
        if (S.id[randomPlot1][individualOne] == randomTaxon2) {
          SppTwoPresentInPlotOne = TRUE;
          sppTwoIndiv = individualOne;
          sppTwoPlot = randomPlot1;
        }
        if (S.id[randomPlot2][individualTwo] == randomTaxon2) {
          SppTwoPresentInPlotTwo = TRUE;
          sppTwoIndiv = individualTwo;
          sppTwoPlot = randomPlot2;
        }
      }
    }

    //if checkboard, swap species between plots
    if (SppOnePresentInPlotOne == SppTwoPresentInPlotTwo) {
      if (SppOnePresentInPlotTwo == SppTwoPresentInPlotOne) {
        if (SppOnePresentInPlotOne != SppTwoPresentInPlotOne) {
          //if it is checkerboard, swap Id and Abund between plots
          //after storing current values and indicate success
          idHolder = S.id[sppOnePlot][sppOneIndiv];
          S.id[sppOnePlot][sppOneIndiv]
            = S.id[sppTwoPlot][sppTwoIndiv];
          S.id[sppTwoPlot][sppTwoIndiv] = idHolder;
          abundHolder = S.abund[sppOnePlot][sppOneIndiv];
          S.abund[sppOnePlot][sppOneIndiv]
            = S.abund[sppTwoPlot][sppTwoIndiv];
          S.abund[sppTwoPlot][sppTwoIndiv] = abundHolder;
        }
      }
    }
  }
}


// --------- Phylogeny Taxa Shuffle ---------------------------
void PhylogenySampleTaxaShuffle(phylo P, sample S, int *attach) {

  //This randomization works by shuffling the terminal taxa labels
  //for species present in the sample
  //and then re-attaching the sample to the phylogeny.
  //Just swaps the Phylo.taxalist labels, so don't need to recalc P.dist
  //Only swaps labels for taxa found in the sample so that species from
  //the phylo are not introduced into the samples.

  int taxon;
  int randomSampleTaxon;
  int foundTaxon;
  int oldAttachTaxon;
  int oldAttachRandomSampleTaxon;
  char tempTaxonName[MAXTAXONLENGTH];

  //swap phylogeny positions of taxa found in the sample
  //for each taxon in the sample...
  //TODO fix sample w/o replacement algorithm - currently broken!
  for (taxon = 0; taxon<S.ntaxa; taxon++) {
    foundTaxon = FALSE;
    while (foundTaxon == FALSE) {
      //...pick another random sample taxon
      randomSampleTaxon = random() % S.ntaxa;
      //can't swap a taxon with itself:
      if (taxon != randomSampleTaxon) {
        //now swap the taxa labels in P.taxalist and P.taxon for those species
        //and swap the attach references as we go so the sample taxa
        //point to the correct phylo taxa.
        oldAttachTaxon = attach[taxon];
        oldAttachRandomSampleTaxon = attach[randomSampleTaxon];

        strcpy(tempTaxonName, P.taxalist[oldAttachTaxon]);
        strcpy(P.taxalist[oldAttachTaxon],
               P.taxalist[oldAttachRandomSampleTaxon]);
        strcpy(P.taxon[P.t2n[oldAttachTaxon]],
               P.taxon[P.t2n[oldAttachRandomSampleTaxon]]);
        strcpy(P.taxalist[oldAttachRandomSampleTaxon], tempTaxonName);
        strcpy(P.taxon[P.t2n[oldAttachRandomSampleTaxon]],
               tempTaxonName);

        //update the attach
        attach[taxon] = oldAttachRandomSampleTaxon;
        attach[randomSampleTaxon] = oldAttachTaxon;

        foundTaxon = TRUE;
      }
    }
  }
}

// --------- Phylogeny Attach Shuffle ---------------------------
void PhylogenyAttachShuffle(phylo P, sample S, int *attach) {

  //This randomization works by shuffling the attach vector connecting
  //sample taxa to their correponding taxa node# on the phylogeny
  //This gives taxa found in the phylo but not the sample a chance to
  //enter the sample, as well as effectively randomizing taxon labels
  //in the phylogeny.

  //attach[sample taxa #] = corresponding phylo taxon node #

  int taxon;
  int tmp;
  int randomPhyloTaxon;
  int *taxaList;

  taxaList = ivector(0, P.termtaxa-1);

  for (taxon = 0; taxon < P.termtaxa; taxon++)
    taxaList[taxon]=taxon;

  //randomize the phylo taxa list
  //NB careful - want terminal taxa(termtaxa?) not internal nodes (ntaxa)
  for (taxon = 0; taxon < P.termtaxa; taxon++) {
    randomPhyloTaxon = random() % P.termtaxa;
    tmp = taxaList[taxon];
    taxaList[taxon] = taxaList[randomPhyloTaxon];
    taxaList[randomPhyloTaxon] = tmp;
  }

  //assign each taxon in sample a new random corresponding phylo taxa node#
  for (taxon = 0; taxon<S.ntaxa; taxon++) {
    attach[taxon] = taxaList[taxon];
    if (Debug)
      printf("%d\t", attach[taxon]);
  }
  if (Debug)
    printf("\n");
}

// ---------- Sample ID + Phylo-Attach Shuffle --------------------------------
void RandomizeSampleTaxaShuffle(sample S) {
  //This method works by assigning each record in the sample file a
  //random sample taxon ID #
  //The sample-phylo attach will need to be subsequently shuffled to
  //give every taxon in the phylo a chance to be included in the samples
  //if the samples don't include all taxa in phylogeny!

  int taxon;
  int sampleNum;
  int recordNum;
  int randomSampleTaxon;
  int tmp;
  int *taxaList;

  taxaList = ivector(0, S.ntaxa-1);

  for (taxon = 0; taxon < S.ntaxa; taxon++)
    taxaList[taxon]=taxon;

  //assign each record in each sample a new random id#
  for (sampleNum = 0; sampleNum<S.nsamples; sampleNum++) {
    if (Debug) printf("sample %d nrecs %d\n",sampleNum,S.srec[sampleNum]);
    if (Debug) {
      for (tmp = 0;tmp<S.srec[sampleNum];tmp++) printf("\trec %d\tid %d\n",tmp,S.id[sampleNum][tmp]);
    }
    for (taxon = 0; taxon < S.srec[sampleNum]; taxon++) {
      randomSampleTaxon = random() % S.ntaxa;
      tmp = taxaList[taxon];
      taxaList[taxon] = taxaList[randomSampleTaxon];
      taxaList[randomSampleTaxon] = tmp;
    }
    for (recordNum = 0; recordNum<S.srec[sampleNum]; recordNum++) {
      if (Debug) printf("samp %d\trec %d\toriginal id %d\t",sampleNum,recordNum,S.id[sampleNum][recordNum]);
      S.id[sampleNum][recordNum] = taxaList[recordNum];
      if (Debug) printf("new id %d\n",S.id[sampleNum][recordNum]);
    }
  }
}

// ---------- OUTPUT SWAPPED MATRIX -------------------------------------------

void OutputSwappedMatrix(phylo InP, sample InS, int SwapMethod) {
  //This function randomizes the phylo or sample and outputs the
  //swapped sample to console in sample (database) format
  //(each line contains plot, stemcount, taxaID)
  //TODO fix so that if phylo not needed, it's not loaded or output
  int plot;
  int taxon;
  int *attach;

  //Attach sample to phylogeny
  attach = ivector(0, InS.ntaxa-1);
  AttachSampleToPhylo(InS, InP, attach);

  if (Debug) {
    for (plot = 0; plot < InS.nsamples; plot++) {
      for (taxon = 0; taxon < InS.srec[plot]; taxon++) {
        printf("%s\t%d\t%s\n", InS.pname[plot], InS.abund[plot][taxon],
               InP.taxalist[attach[InS.id[plot][taxon]]]);
      }
    }
  }

  //Swap the matrix or phylogeny
  switch (SwapMethod) {
  case 0:
    //SwapMethod 0 = shuffle taxa labels on phylogeny
    //(sample remains unshuffled)
    PhylogenyAttachShuffle(InP, InS, attach);
    break;
  case 1:
    //SwapMethod 1 = samples become random draws from sample taxa list
    //(maintains sample species richnesses but not species frequencies)
    //(Species present in the phylogeny but not the sample will not
    //be represented in the shuffled samples)
    RandomizeSampleTaxaShuffle(InS);
    break;
  case 2:
    //SwapMethod 2 = samples become random draws from phylogeny
    //(maintains sample species richnesses but not species frequencies)
    //(also shuffles link between phylo and sample taxa to give
    //species present in phylo but not sample chance to be represented)
    //Note limitation of this way of shuffling is that if phylogeny
    //contains more taxa than sample, each shuffle will only
    //fill the samples with sample # taxa, not phylo # taxa
    RandomizeSampleTaxaShuffle(InS);
    PhylogenyAttachShuffle(InP, InS, attach);
    break;
  case 3:
    //SwapMethod 3 = independent checkerboard swap of sample matrix
    //(maintains species frequencies and sample species richnesses)
    IndependentSwap(InS, SWAPS);
    break;
    //          case x:
    //              //SwapMethod x = independent swap & phylogeny shuffle
    //              //(combine swapmethods 0 & 1)
    //              GotelliSwap(InS,SWAPS);
    //              PhylogenyAttachShuffle(InP,InS,attach);
    //              break;
    //          case y:
    //              //SwapMethod y = shuffle phylo taxa labels and reattach phylogeny
    //              //to sample
    //              //This is similar to swapmethod 0 (shuffle phylogeny attach), with
    //              //the difference being that by shuffling labels and then re-attaching
    //              //sample to phylo, taxa in phylo but not sample will not be included
    //              //in randomized samples (whereas they will be included using method 0)
    //              PhylogenyTaxaShuffle(InP,InS,attach);
    //              break;
  default:
    printf("Please use -m command line switch to specify a randomization method.\n");
    printf("See documentation for a list of possible null models.\n");
    exit(EXIT_FAILURE);
    break;
  }

  //Output the swapped matrix to console
  //Note - using attach to insert taxa names from phylogeny in case
  //sample-phylo attach has been shuffled.
  for (plot = 0; plot < InS.nsamples; plot++) {
    for (taxon = 0; taxon < InS.srec[plot]; taxon++) {
      printf("%s\t%d\t%s\n", InS.pname[plot], InS.abund[plot][taxon],
             InP.taxalist[attach[InS.id[plot][taxon]]]);
    }
  }
}

// --------- MeanDistance -----------------------------------------------------

double MeanDistance(phylo InP, sample InS, int *attach, int plot,
                    int abundWeighted) {

  if (abundWeighted == 1) {
    //equals sum (abundi * abundj * distij) / sum (abundi * abundj)
    //conspecifics are counted as part of the average distance (we count individuals here)

    int j, k;
    double distSum = 0.0;
    long abundSum = 0;

    for (j = 0; j < InS.srec[plot]; j++) {
      for (k = 0; k < InS.srec[plot]; k++) {
        //Sum distances among all pairs of individuals in the plot
        distSum
          += InS.abund[plot][j] * InS.abund[plot][k]
          * InP.dist[InP.t2n[attach[InS.id[plot][j]]]][InP.t2n[attach[InS.id[plot][k]]]];
        abundSum += InS.abund[plot][j] * InS.abund[plot][k];
      }
    }
    //calculate weighted mean pairwise distance among taxa for this plot
    return distSum / (double) abundSum;
  } else {
    //Non-abundance weighted
    //For j != k
    int j, k;
    double plotSum = 0.0;
    long pairCount = 0;

    for (j = 0; j < InS.srec[plot]-1; j++) {
      for (k = j+1; k < InS.srec[plot]; k++) {
        //Sum distances among all pairs of taxa in the plot
        plotSum
          += InP.dist[InP.t2n[attach[InS.id[plot][j]]]][InP.t2n[attach[InS.id[plot][k]]]];
        pairCount++;
      }
    }
    //calculate mean pairwise distance among taxa for this plot
    if (Debug)
      printf("%f\t%ld\n", plotSum, pairCount);
    return plotSum/(double)pairCount;
  }
}

// ---------- MeanMinimumDistance ---------------------------------------------

double MeanMinimumDistance(phylo InP, sample InS, int *attach, int plot,
                           int abundWeighted) {
  int j, k;
  double thisDist = 0.0;
  double minDist = 0.0;
  double sumMinDist = 0.0;
  double abundSum = 0.0;

  if (InS.srec[plot]==1) {
    return 0.0;
  }
  else {
    for (j = 0; j < InS.srec[plot]; j++) {
      minDist = 99999999.9; //needs to be higher than expected min distance

      for (k = 0; k < InS.srec[plot]; k++) {
        if (k != j) {
          //find lowest pairwise distance between taxa
          thisDist
            = InP.dist[InP.t2n[attach[InS.id[plot][j]]]][InP.t2n[attach[InS.id[plot][k]]]];
          if (thisDist < minDist)
            minDist = thisDist;
        }
      }
      if (abundWeighted == 1) {
        minDist = (double)InS.abund[plot][j] * minDist;
        abundSum += (double)InS.abund[plot][j];
      }
      sumMinDist += minDist;
    }

    //return the average minimum pairwise distance among taxa in the plot
    if (abundWeighted == 1) {
      return sumMinDist/(double)abundSum;
    } else {
      return sumMinDist/(double)InS.srec[plot];
    }
  }
}

// Community distance with null comparison (abundance weighted)
void CommunityDistanceNull(phylo InP, sample InS, int SwapMethod, int abundWeighted, int nullTest) {
  int j, k, plot1, plot2, runNum;
  int *attach;
  float **comDist;
  float ***comDistRnd;
  float **comDistRndMean;
  float **comDistRndSD;
  comDist = matrix(0, InS.nsamples-1, 0, InS.nsamples-1);
  comDistRnd = f3tensor(0, 1, 0, 1, 0, 1);
  comDistRndMean = matrix(0, 1, 0, 1);
  comDistRndSD = matrix(0, 1, 0, 1);

  if (nullTest == 1) {

    comDistRnd = f3tensor(0, RUNS-1, 0, InS.nsamples-1, 0, InS.nsamples-1);
    comDistRndMean = matrix(0, InS.nsamples-1, 0, InS.nsamples-1);
    comDistRndSD = matrix(0, InS.nsamples-1, 0, InS.nsamples-1);

    //initialize matrices
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        comDistRndMean[plot1][plot2] = 0.0;
        comDistRndSD[plot1][plot2] = 0.0;
      }
    }

  }

  //Create and fill matrix of distances among all nodes in the phylogeny
  InP.dist = matrix(0, InP.nnodes-1, 0, InP.nnodes-1);
  DistMatrix(InP);

  //Attach sample to phylogeny
  attach = ivector(0, InS.ntaxa-1);
  AttachSampleToPhylo(InS, InP, attach);

  //Calculate observed comdist
  if (abundWeighted == 1) {
    //equals sum (pabundi * pabundj * distij) / sum (pabundi * pabundj)
    float distSum = 0.0;
    float abundSum = 0.0;
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        distSum = 0.0;
        abundSum = 0.0;
        for (j = 0; j < InS.srec[plot1]; j++) {
          for (k = 0; k < InS.srec[plot2]; k++) {
            //Sum distances among all pairs of taxa in the 2 plots
            distSum += InS.pabund[plot1][j] * InS.pabund[plot2][k]
              * TaxaDist(InP, InS.id[plot1][j],
                         InS.id[plot2][k], attach);
            abundSum += InS.pabund[plot1][j] * InS.pabund[plot2][k];
          }
        }
        //calculate weighted mean pairwise distance among taxa
        comDist[plot1][plot2] = distSum / abundSum;
      }
    }
  } else {
    //Non-abundance weighted
    //NB Within-community distance calculated differently than among-community distance
    //Within-community distance = MPD (does not include conspecifics in calculation)
    //Among-community distance does include conspecifics in calculation of mean distance
    float plotSum = 0.0;
    long pairCount = 0;
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        plotSum = 0.0;
        pairCount = 0;
        for (j = 0; j < InS.srec[plot1]; j++) {
          for (k = 0; k < InS.srec[plot2]; k++) {
            //Sum distances among all pairs of taxa in the 2 plots
            plotSum += TaxaDist(InP, InS.id[plot1][j],
                                InS.id[plot2][k], attach);
            pairCount++;
          }
        }
        //calculate weighted mean pairwise distance among taxa
        comDist[plot1][plot2] = plotSum / (float) pairCount;
        //if (Debug == 1) printf("plot1\t%d\tplot2\t%d\tcomdist\t%f\n",plot1,plot2,comDist[plot1][plot2]);

      }
      //Within-community distance = MPD
      comDist[plot1][plot1] = MeanDistance(InP, InS, attach, plot1, 0);
    }
  }

  if (nullTest == 1)
    {
      // if BURNIN > 0 and independent or trial swap null, randomize and discard (burn in)
      if ((BURNIN > 0) && (SwapMethod==3 || SwapMethod==4)) {
        //burn in using appropriate algorithm
        if (SwapMethod == 3) {
          IndependentSwap(InS, BURNIN);
        } else {
          TrialSwap(InS, BURNIN);
        }
      }

      //Calculate comdists for randomized data
      for (runNum = 0; runNum < RUNS; runNum++) {
        //initialize matrices
        for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
          for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
            comDistRnd[runNum][plot1][plot2] = 0.0;
          }
        }

        //Swap the matrix or phylogeny
        switch (SwapMethod) {
        case 0:
          //SwapMethod 0 = shuffle taxa labels on phylogeny
          //(sample remains unshuffled)
          PhylogenyAttachShuffle(InP, InS, attach);
          break;
        case 1:
          //SwapMethod 1 = samples become random draws from sample taxa list
          //(maintains sample species richnesses but not species frequencies)
          //(Species present in the phylogeny but not the sample will not
          //be represented in the shuffled samples)
          RandomizeSampleTaxaShuffle(InS);
          break;
        case 2:
          //SwapMethod 2 = samples become random draws from phylogeny
          //(maintains sample species richnesses but not species frequencies)
          //(also shuffles link between phylo and sample taxa to give
          //species present in phylo but not sample chance to be represented)
          //Note limitation of this way of shuffling is that if phylogeny
          //contains more taxa than sample, each shuffle will only
          //fill the sample with sample # taxa, not phylo # taxa
          RandomizeSampleTaxaShuffle(InS);
          PhylogenyAttachShuffle(InP, InS, attach);
          break;
        case 3:
          //SwapMethod 3 = independent checkerboard swap of sample matrix
          //(maintains species frequencies and sample species richnesses)
          IndependentSwap(InS, SWAPS);
          break;
        case 4:
          //SwapMethod 4 = trial swap of sample matrix
          //(maintains species frequencies and sample species richnesses)
          TrialSwap(InS, SWAPS);
          break;
        default:
          printf("Please use -m command line switch to specify a randomization method.\n");
          printf("See documentation for a list of possible null models.\n");
          exit(EXIT_FAILURE);
          break;
        }

        if (abundWeighted == 1) {
          //equals sum (pabundi * pabundj * distij) / sum (pabundi * pabundj)
          double distSum = 0.0;
          double abundSum = 0.0;
          for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
            for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
              distSum = 0.0;
              abundSum = 0.0;
              for (j = 0; j < InS.srec[plot1]; j++) {
                for (k = 0; k < InS.srec[plot2]; k++) {
                  //Sum distances among all pairs of taxa in the 2 plots
                  distSum += InS.pabund[plot1][j] * InS.pabund[plot2][k]
                    * TaxaDist(InP, InS.id[plot1][j],
                               InS.id[plot2][k], attach);
                  abundSum += InS.pabund[plot1][j] * InS.pabund[plot2][k];
                }
              }
              //calculate weighted mean pairwise distance among taxa
              comDistRnd[runNum][plot1][plot2] = distSum / abundSum;
              comDistRndMean[plot1][plot2] += distSum / abundSum;
            }
          }
        } else {
          //Non-abundance weighted
          //NB Within-community distance calculated differently than among-community distance
          //Within-community distance = MPD (does not include conspecifics in calculation)
          //Among-community distance does include conspecifics in calculation of mean distance
          //TODO fixme incorrect weighting of plots with unequal abundance?
          double plotSum = 0.0;
          long pairCount = 0;
          for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
            for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
              plotSum = 0.0;
              pairCount = 0;
              for (j = 0; j < InS.srec[plot1]; j++) {
                for (k = 0; k < InS.srec[plot2]; k++) {
                  //Sum distances among all pairs of taxa in the 2 plots
                  plotSum += TaxaDist(InP, InS.id[plot1][j],
                                      InS.id[plot2][k], attach);
                  pairCount++;
                }
              }
              //calculate weighted mean pairwise distance among taxa
              if (plot1 != plot2) {
                comDistRnd[runNum][plot1][plot2] = plotSum / (float) pairCount;
                comDistRndMean[plot1][plot2] += plotSum / (float) pairCount;
              }

            }
            //Within-community distance = MPD
            comDistRnd[runNum][plot1][plot1] = MeanDistance(InP, InS, attach, plot1, 0);
            comDistRndMean[plot1][plot1] += MeanDistance(InP, InS, attach, plot1, 0);
          }
        }
      }

      //calculate randomization results - mean, SD of random comDists
      for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
        for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
          comDistRndMean[plot1][plot2] /= (double)RUNS;
          for (runNum = 0; runNum < RUNS; runNum++) {
            //sum squares
            comDistRndSD[plot1][plot2] += pow(comDistRndMean[plot1][plot2] - comDistRnd[runNum][plot1][plot2], 2)
              /(double)(RUNS - 1);
          }
          comDistRndSD[plot1][plot2] = sqrt(comDistRndSD[plot1][plot2]);
        }
      }

    }

  //summarize results
  if (nullTest == 1) {
    printf("#ObservedComdist\n");
  }
  printf(".");
  for (plot1 = 0; plot1 < InS.nsamples; plot1++)
    printf("\t%s", InS.pname[plot1]);
  printf("\n");
  for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
    printf("%s", InS.pname[plot1]);
    for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
      printf("\t%f", comDist[plot1][plot2]);
    }
    printf("\n");
  }

  if (nullTest == 1) {
    printf("\n#RandomMeanComdist\n");
    printf(".");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++)
      printf("\t%s", InS.pname[plot1]);
    printf("\n");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      printf("%s", InS.pname[plot1]);
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        printf("\t%f", comDistRndMean[plot1][plot2]);
      }
      printf("\n");
    }

    printf("\n#RandomSDComdist\n");
    printf(".");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++)
      printf("\t%s", InS.pname[plot1]);
    printf("\n");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      printf("%s", InS.pname[plot1]);
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        printf("\t%f", comDistRndSD[plot1][plot2]);
      }
      printf("\n");
    }

    printf("\n#NRIComdist\n");
    printf(".");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++)
      printf("\t%s", InS.pname[plot1]);
    printf("\n");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      printf("%s", InS.pname[plot1]);
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        printf("\t%f", (comDistRndMean[plot1][plot2] - comDist[plot1][plot2])/comDistRndSD[plot1][plot2] );
      }
      printf("\n");
    }

    free_matrix(comDistRndMean, 0, InS.nsamples-1, 0, InS.nsamples-1);
    free_matrix(comDistRndSD, 0, InS.nsamples-1, 0, InS.nsamples-1);
    free_f3tensor(comDistRnd, 0, RUNS-1, 0, InS.nsamples-1, 0, InS.nsamples-1);

  }

  //free memory
  free_matrix(comDist, 0, InS.nsamples-1, 0, InS.nsamples-1);
  free_ivector(attach, 0, InS.ntaxa-1);
}


// Community distance with null comparison (abundance weighted)
void CommunityDistanceNNNull(phylo InP, sample InS, int SwapMethod, int abundWeighted, int nullTest) {

  float thisDist = 0.0;
  float minDist = 0.0;
  float abundSum = 0.0;
  int j, k, plot1, plot2, runNum;
  int *attach;
  float **comDist;
  float ***comDistRnd;
  float **comDistRndMean;
  float **comDistRndSD;

  comDist = matrix(0, InS.nsamples-1, 0, InS.nsamples-1);
  comDistRnd = f3tensor(0, 1, 0, 1, 0, 1);
  comDistRndMean = matrix(0, 1, 0, 1);
  comDistRndSD = matrix(0, 1, 0, 1);

  if (nullTest == 1) {
    comDistRnd = f3tensor(0, RUNS-1, 0, InS.nsamples-1, 0, InS.nsamples-1);
    comDistRndMean = matrix(0, InS.nsamples-1, 0, InS.nsamples-1);
    comDistRndSD = matrix(0, InS.nsamples-1, 0, InS.nsamples-1);

    //initialize matrices
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        comDistRndMean[plot1][plot2] = 0.0;
        comDistRndSD[plot1][plot2] = 0.0;
      }
    }
  }

  //Create and fill matrix of distances among all nodes in the phylogeny
  InP.dist = matrix(0, InP.nnodes-1, 0, InP.nnodes-1);
  DistMatrix(InP);

  //Attach sample to phylogeny
  attach = ivector(0, InS.ntaxa-1);
  AttachSampleToPhylo(InS, InP, attach);

  for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
    for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
      if (plot1 != plot2) {
        comDist[plot1][plot2] = 0.0;
        abundSum = 0.0;
        for (j = 0; j < InS.srec[plot1]; j++) {
          minDist = 99999999.9; //needs to be higher than expected min distance
          for (k = 0; k < InS.srec[plot2]; k++) {
            //find lowest pairwise distance between taxa
            thisDist = InP.dist[InP.t2n[attach[InS.id[plot1][j]]]][InP.t2n[attach[InS.id[plot2][k]]]];
            if (thisDist < minDist)
              minDist = thisDist;
          }
          if (abundWeighted == 1) {
            minDist = (float)InS.pabund[plot1][j] * minDist;
            abundSum += (float)InS.pabund[plot1][j];
          }
          else
            {
              minDist = (1.0 / (double)InS.srec[plot1]) * minDist;
              abundSum += (1.0 / (double)InS.srec[plot1]);
            }
          comDist[plot1][plot2] += minDist;
        }
        //XXX comment the following out for non-symmetric distances
        for (j = 0; j < InS.srec[plot2]; j++) {
          minDist = 99999999.9; //needs to be higher than expected min distance
          for (k = 0; k < InS.srec[plot1]; k++) {
            //find lowest pairwise distance between taxa
            thisDist = InP.dist[InP.t2n[attach[InS.id[plot2][j]]]][InP.t2n[attach[InS.id[plot1][k]]]];
            if (thisDist < minDist)
              minDist = thisDist;
          }
          if (abundWeighted == 1) {
            minDist = (float)InS.pabund[plot2][j] * minDist;
            abundSum += (float)InS.pabund[plot2][j];
          }
          else
            {
              minDist = (1.0 / (float)InS.srec[plot2]) * minDist;
              abundSum += (1.0 / (float)InS.srec[plot2]);
            }
          comDist[plot1][plot2] += minDist;
        }
        //XXX
        //calculate weighted mean pairwise distance among taxa
        comDist[plot1][plot2] = comDist[plot1][plot2]/(float)abundSum;
      }
      else comDist[plot1][plot2] = MeanMinimumDistance(InP, InS, attach, plot1, abundWeighted);
    }
  }

  if (nullTest == 1) {

    // if BURNIN > 0 and independent or trial swap null, randomize and discard (burn in)
    if ((BURNIN > 0) && (SwapMethod==3 || SwapMethod==4)) {
      //burn in using appropriate algorithm
      if (SwapMethod == 3) {
        IndependentSwap(InS, BURNIN);
      } else {
        TrialSwap(InS, BURNIN);
      }
    }

    //Calculate comdists for randomized data
    for (runNum = 0; runNum < RUNS; runNum++) {

      //Swap the matrix or phylogeny
      switch (SwapMethod) {
      case 0:
        //SwapMethod 0 = shuffle taxa labels on phylogeny
        //(sample remains unshuffled)
        PhylogenyAttachShuffle(InP, InS, attach);
        break;
      case 1:
        //SwapMethod 1 = samples become random draws from sample taxa list
        //(maintains sample species richnesses but not species frequencies)
        //(Species present in the phylogeny but not the sample will not
        //be represented in the shuffled samples)
        RandomizeSampleTaxaShuffle(InS);
        break;
      case 2:
        //SwapMethod 2 = samples become random draws from phylogeny
        //(maintains sample species richnesses but not species frequencies)
        //(also shuffles link between phylo and sample taxa to give
        //species present in phylo but not sample chance to be represented)
        //Note limitation of this way of shuffling is that if phylogeny
        //contains more taxa than sample, each shuffle will only
        //fill the sample with sample # taxa, not phylo # taxa
        RandomizeSampleTaxaShuffle(InS);
        PhylogenyAttachShuffle(InP, InS, attach);
        break;
      case 3:
        //SwapMethod 3 = independent checkerboard swap of sample matrix
        //(maintains species frequencies and sample species richnesses)
        IndependentSwap(InS, SWAPS);
        break;
      case 4:
        //SwapMethod 4 = trial swap of sample matrix
        //(maintains species frequencies and sample species richnesses)
        TrialSwap(InS, SWAPS);
        break;
      default:
        printf("Please use -m command line switch to specify a randomization method.\n");
        printf("See documentation for a list of possible null models.\n");
        exit(EXIT_FAILURE);
        break;
      }

      for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
        for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
          if (plot1 != plot2) {
            comDistRnd[runNum][plot1][plot2] = 0.0;
            abundSum = 0.0;
            for (j = 0; j < InS.srec[plot1]; j++) {
              minDist = 99999999.9; //needs to be higher than expected min distance
              for (k = 0; k < InS.srec[plot2]; k++) {
                //find lowest pairwise distance between taxa
                thisDist = InP.dist[InP.t2n[attach[InS.id[plot1][j]]]][InP.t2n[attach[InS.id[plot2][k]]]];
                if (thisDist < minDist)
                  minDist = thisDist;
              }
              if (abundWeighted == 1) {
                minDist = (float)InS.abund[plot1][j] * minDist;
                abundSum += (float)InS.abund[plot1][j];
              }
              else
                {
                  minDist = (1.0 / (float)InS.srec[plot1]) * minDist;
                  abundSum += (1.0 / (double)InS.srec[plot1]);
                }
              comDistRnd[runNum][plot1][plot2] += minDist;
            }
            //XXX comment the following out for non-symmetric distances
            for (j = 0; j < InS.srec[plot2]; j++) {
              minDist = 99999999.9; //needs to be higher than expected min distance
              for (k = 0; k < InS.srec[plot1]; k++) {
                //find lowest pairwise distance between taxa
                thisDist = InP.dist[InP.t2n[attach[InS.id[plot2][j]]]][InP.t2n[attach[InS.id[plot1][k]]]];
                if (thisDist < minDist)
                  minDist = thisDist;
              }
              if (abundWeighted == 1) {
                minDist = (double)InS.abund[plot2][j] * minDist;
                abundSum += (double)InS.abund[plot2][j];
              }
              else
                {
                  minDist = (1.0 / (double)InS.srec[plot2]) * minDist;
                  abundSum += (1.0 / (double)InS.srec[plot2]);
                }
              comDistRnd[runNum][plot1][plot2] += minDist;
            }
            //XXX
            //calculate weighted mean pairwise distance among taxa
            comDistRnd[runNum][plot1][plot2] = comDistRnd[runNum][plot1][plot2]/(double)abundSum;
            comDistRndMean[plot1][plot2] += comDistRnd[runNum][plot1][plot2];
          } else {
            comDistRnd[runNum][plot1][plot2] = MeanMinimumDistance(InP, InS, attach, plot1, abundWeighted);
            comDistRndMean[plot1][plot2] += comDistRnd[runNum][plot1][plot2];
          }
        }
      }
    }

    //calculate randomization results - mean, SD of random comDists
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        comDistRndMean[plot1][plot2] /= (double)RUNS;
        for (runNum = 0; runNum < RUNS; runNum++) {
          //sum squares
          comDistRndSD[plot1][plot2] += pow(comDistRndMean[plot1][plot2] - comDistRnd[runNum][plot1][plot2], 2)
            /(double)(RUNS - 1);
          //if (Debug) printf("%f\t",comDistRnd[runNum][plot1][plot2]);
        }
        comDistRndSD[plot1][plot2] = sqrt(comDistRndSD[plot1][plot2]);
      }
    }
  }

  //summarize results
  if (nullTest == 1) {
    printf("#ObservedComdistNT\n");
  }
  printf(".");
  for (plot1 = 0; plot1 < InS.nsamples; plot1++)
    printf("\t%s", InS.pname[plot1]);
  printf("\n");
  for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
    printf("%s", InS.pname[plot1]);
    for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
      printf("\t%f", comDist[plot1][plot2]);
    }
    printf("\n");
  }

  if (nullTest == 1) {
    printf("\n#RandomMeanComdistNT\n");
    printf(".");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++)
      printf("\t%s", InS.pname[plot1]);
    printf("\n");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      printf("%s", InS.pname[plot1]);
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        printf("\t%f", comDistRndMean[plot1][plot2]);
      }
      printf("\n");
    }

    printf("\n#RandomSDComdistNT\n");
    printf(".");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++)
      printf("\t%s", InS.pname[plot1]);
    printf("\n");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      printf("%s", InS.pname[plot1]);
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        printf("\t%f", comDistRndSD[plot1][plot2]);
      }
      printf("\n");
    }

    printf("\n#NTIComdistNT\n");
    printf(".");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++)
      printf("\t%s", InS.pname[plot1]);
    printf("\n");
    for (plot1 = 0; plot1 < InS.nsamples; plot1++) {
      printf("%s", InS.pname[plot1]);
      for (plot2 = 0; plot2 < InS.nsamples; plot2++) {
        printf("\t%f", (comDistRndMean[plot1][plot2] - comDist[plot1][plot2])/comDistRndSD[plot1][plot2] );
      }
      printf("\n");
    }

    free_matrix(comDistRndMean, 0, InS.nsamples-1, 0, InS.nsamples-1);
    free_matrix(comDistRndSD, 0, InS.nsamples-1, 0, InS.nsamples-1);
    free_f3tensor(comDistRnd, 0, RUNS-1, 0, InS.nsamples-1, 0, InS.nsamples-1);
  }

  //free memory
  free_matrix(comDist, 0, InS.nsamples-1, 0, InS.nsamples-1);
  free_ivector(attach, 0, InS.ntaxa-1);

}


// ---------- PhyloDiversity ---------------------------------------------
void PhyloDiversity(phylo P, sample S) {
  //Calculates community phylogenetic diversity
  //Notation/formulae follow Hardy & Senterre 2007
  //Currently in flux - don't use!

  int i, j, k, l;
  int plot;
  double sumDist;
  long countDist;
  double sampDist;
  long sampCount;
  int *attach;

  double N = (double)S.nsamples;

  //diversity coefficients
  double *Di;
  double *Dp;
  double *DeltaP;
  double DiS, DpS, DeltaS;
  double DiT, DpT, DeltaT;
  double Ist, Pst;
  double PIst;
  //unbiased estimators of diversity coefficients
  double *DiHat;
  double *DpHat;
  double DiSHat, DpSHat;
  double DiTHat, DpTHat, DeltaTHat;
  double IstHat, PstHat;
  double PIstHat;

  //allocate vectors and matrices
  Di = dvector(0, S.nsamples-1);
  Dp = dvector(0, S.nsamples-1);
  DiHat = dvector(0, S.nsamples-1);
  DpHat = dvector(0, S.nsamples-1);
  DeltaP = dvector(0, S.nsamples-1);

  //Attach sample to phylogeny
  attach = ivector(0, S.ntaxa-1);
  AttachSampleToPhylo(S, P, attach);

  //Create and fill matrix of distances among all nodes in the phylogeny
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  //calculate within-site diversity
  DiS = 0.0;
  DpS = 0.0;
  DiSHat = 0.0;
  DpSHat = 0.0;
  DeltaS = 0.0;
  for (plot = 0; plot < S.nsamples; plot++) {
    Di[plot] = 0.0;
    DiHat[plot] = 0.0;
    Dp[plot] = 0.0;
    DpHat[plot] = 0.0;
    DeltaP[plot] = 0.0;
    countDist = 0;
    for (i = 0; i < S.srec[plot]; i++) {
      for (j = 0; j < S.srec[plot]; j++) {
        //if (i != j) {
        Dp[plot] += (TaxaDist(P, S.id[plot][i], S.id[plot][j], attach)
                     / 2.0) * S.pabund[plot][i]* S.pabund[plot][j];
        if (i != j) {
          Di[plot] += S.pabund[plot][i] * S.pabund[plot][j];
          DeltaP[plot] += (TaxaDist(P, S.id[plot][i], S.id[plot][j],
                                    attach) / 2.0);
          countDist++;
        }
        //}
      }
    }
    DiHat[plot] = Di[plot] * ( (double)S.irec[plot] / ((double)S.irec[plot]
                                                       - 1.0));
    DpHat[plot] = Dp[plot] * ( (double)S.irec[plot] / ((double)S.irec[plot]
                                                       - 1.0));
    DeltaP[plot] /= (double)countDist;
  }

  //calculate average within-site diversities weighted by site abund/richness
  double sumWt = 0;
  double sumWtDelta = 0;
  for (k = 0; k < S.nsamples; k++) {
    double n = (double)S.irec[k];
    sumWt += (n * (n-1));
    DiS += Di[k];
    DpS += Dp[k];
    DiSHat += DiHat[k] * n * (n-1);
    DpSHat += DpHat[k] * n * (n-1);
    //fix:
    n = (double)S.srec[k];
    sumWtDelta += (n * (n-1));
    DeltaS += DeltaP[k] * (n * (n-1));
  }
  DiS /= S.nsamples;
  DpS /= S.nsamples;
  DiSHat /= sumWt;
  DpSHat /= sumWt;
  DeltaS /= sumWtDelta;

  double DiW = DiSHat;
  double DpW = DpSHat;
  double DeltaW = DeltaS;

  //calculate among-site diversity
  double **DiIJ;
  double **DpIJ;
  double **DeltaPIJ;
  DiIJ = dmatrix(0, S.nsamples-1, 0, S.nsamples-1);
  DpIJ = dmatrix(0, S.nsamples-1, 0, S.nsamples-1);
  DeltaPIJ = dmatrix(0, S.nsamples-1, 0, S.nsamples-1);
  unsigned long countDiff;
  for (k = 0; k < S.nsamples; k++) {
    for (l = 0; l < S.nsamples; l++) {
      DiIJ[k][l] = 0.0;
      DpIJ[k][l] = 0.0;
      DeltaPIJ[k][l] = 0.0;
      countDiff = 0;
      if (k != l) {
        for (i = 0; i < S.srec[k]; i++) {
          for (j = 0; j < S.srec[l]; j++) {
            if (S.id[k][i] != S.id[l][j]) {
              DiIJ[k][l] += S.pabund[k][i] * S.pabund[l][j];
              DpIJ[k][l] += S.pabund[k][i] * S.pabund[l][j]
                * (TaxaDist(P, S.id[k][i], S.id[l][j],
                            attach) / 2.0 );
              DeltaPIJ[k][l] += (TaxaDist(P, S.id[k][i],
                                          S.id[l][j], attach) / 2.0 );
              countDiff++;
            }
          }
        }
      }
      DeltaPIJ[k][l] /= (double)countDiff;
    }
  }

  //calculate average among-site diversities weighted by site abund/richness
  double DiA = 0.0;
  double DpA = 0.0;
  double DeltaA = 0.0;
  sumWt = 0;
  sumWtDelta = 0;
  for (k = 0; k < S.nsamples; k++) {
    for (l = 0; l < S.nsamples; l++) {
      if (k != l) {
        double n1 = (double)S.irec[k];
        double n2 = (double)S.irec[l];
        double N = (double)S.totabund;
        sumWt += (n1 * (N-n1)) * (n2 * (N-n2));
        DiA += DiIJ[k][l] * (n1 * (N-n1)) * (n2 * (N-n2));
        DpA += DpIJ[k][l] * (n1 * (N-n1)) * (n2 * (N-n2));
        //fix:
        n1 = (double)S.srec[k];
        n2 = (double)S.srec[l];
        N = (double)S.ntaxa;
        sumWtDelta += (n1 * (N-n1)) * (n2 * (N-n2));
        DeltaA += DeltaPIJ[k][l] * (n1 * (N-n1)) * (n2 * (N-n2));
      }
    }
  }
  DiA /= sumWt;
  DpA /= sumWt;
  DeltaA /= sumWtDelta;

  //calculate total diversities
  DiT = 0.0;
  DpT = 0.0;
  sumDist = 0.0;
  countDist = 0.0;
  for (i = 0; i < S.ntaxa; i++) {
    for (j = 0; j < S.ntaxa; j++) {
      if (i != j) {
        DiT += S.psppabund[i] * S.psppabund[j];
        DpT += (TaxaDist(P, i, j, attach) / 2.0) * S.psppabund[i]
          * S.psppabund[j];
        sumDist += (TaxaDist(P, i, j, attach) / 2.0);
        countDist++;
      }
    }
  }
  DeltaT = sumDist / (double)countDist;
  //    Ist = (DiT - DiS) / DiT;
  //    Pst = (DpT - DpS) / DpT;
  //    PIst = (DeltaT - DeltaS) / DeltaT;

  Ist = 1 - (DiW/DiA);
  Pst = 1 - (DpW/DpA);
  PIst = 1 - (DeltaW/DeltaA);

  //calculate unbiased estimators of total diversity
  DiTHat = 0.0;
  DpTHat = 0.0;
  sumDist = 0.0;
  countDist = 0.0;
  for (i = 0; i < S.ntaxa; i++) {
    for (j = 0; j < S.ntaxa; j++) {
      if (i != j) {
        DiTHat += S.psppabund[i] * S.psppabund[j];
        DpTHat += (TaxaDist(P, i, j, attach) / 2.0) * S.psppabund[i]
          * S.psppabund[j];
      }
    }
  }
  IstHat = (DiTHat - DiSHat) / DiTHat;
  PstHat = (DpTHat - DpSHat) / DpTHat;
  //calculate DeltaTHat
  for (k = 0; k < S.nsamples; k++) {
    for (l = 0; l < S.nsamples; l++) {
      if (k != l) {
        sampDist = 0.0;
        sampCount = 0;
        for (i = 0; i < S.srec[k]; i++) {
          for (j = 0; j < S.srec[l]; j++) {
            if (i != j) {
              sampDist += (TaxaDist(P, S.id[k][i], S.id[l][j],
                                    attach) / 2.0);
              sampCount++;
            }
          }
        }
        sumDist += (sampDist / (double)sampCount);
      }
    }
  }
  DeltaTHat = (1.0 / (N * (N - 1.0))) * sumDist;
  PIstHat = (DeltaTHat - DeltaS) / DeltaTHat;

  //write results
  printf("#----- SUMMARY\n");
  printf("DiT\t%f\tDiTHat\t%f\n", DiT, DiTHat);
  printf("DpT\t%f\tDpTHat\t%f\n", DpT, DpTHat);
  printf("DiS\t%f\tDiSHat\t%f\n", DiS, DiSHat);
  printf("DpS\t%f\tDpSHat\t%f\n", DpS, DpSHat);
  printf("Ist\t%f\tIstHat\t%f\n", Ist, IstHat);
  printf("Pst\t%f\tPstHat\t%f\n", Pst, PstHat);
  printf("Pst-Ist\t%f\tPstHat-IstHat\t%f\n", Pst-Ist, PstHat-IstHat);
  printf("DeltaT\t%f\tDeltaTHat\t%f\n", DeltaT, DeltaTHat);
  printf("DeltaS\t%f\n", DeltaS);
  printf("PIst\t%f\tPIstHat\t%f\n", PIst, PIstHat);
  printf("\n");
  //printf("DiA\t%f\tDpA\t%f\tDeltaA\t%f\n",DiA,DpA,DeltaA);
  printf("DiW\t%f\tDiA\t%f\tIst\t%f\n", DiW, DiA, Ist);
  printf("DpW\t%f\tDpA\t%f\tPst\t%f\n", DpW, DpA, Pst);
  printf("DeltaW\t%f\tDeltaA\t%f\tPIst\t%f\n", DeltaW, DeltaA, PIst);

  printf("#----- WITHIN\n");
  printf("Plot\tNSpp\tNIndiv\tDi\tDiHat\tDp\tDpHat\tDeltaP\n");
  for (plot = 0; plot < S.nsamples; plot++) {
    printf("%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n", S.pname[plot], S.srec[plot],
           S.irec[plot], Di[plot], DiHat[plot], Dp[plot], DpHat[plot],
           DeltaP[plot]);
  }
  printf("#----- AMONG - Di.kl\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      printf("\t%f", DiIJ[k][l]);
    }
    printf("\n");
  }
  printf("#----- AMONG - Dp.kl\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      printf("\t%f", DpIJ[k][l]);
    }
    printf("\n");
  }
  printf("#----- AMONG - DELTAp.kl\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      printf("\t%f", DeltaPIJ[k][l]);
    }
    printf("\n");
  }

  //TODO free vectors and matrices from memory

}

// ---------- RaoDiversity ---------------------------------------------
void RaoDiversity(phylo P, sample S) {
  //Calculates Rao's quadratic entropy
  //(alpha/beta diversity weighted/unweighted by phylogenetic distances)
  //Notation follows Rao 1982, Champeley and Chessel 2004
  int i, j, k, l;
  int plot;
  double sumDist;
  long countDist;
  int *attach;

  //diversity coefficients
  double *D;
  double *Dp;
  double **H;
  double **Hp;
  double Dalpha, Dbeta, Dtotal;
  double Dalphap, Dbetap, Dtotalp;
  double Fst, Fstp;

  //allocate vectors and matrices
  D = dvector(0, S.nsamples-1);
  Dp = dvector(0, S.nsamples-1);
  H = dmatrix(0, S.nsamples-1, 0, S.nsamples-1);
  Hp = dmatrix(0, S.nsamples-1, 0, S.nsamples-1);

  //Attach sample to phylogeny
  attach = ivector(0, S.ntaxa-1);
  AttachSampleToPhylo(S, P, attach);

  //Create and fill matrix of distances among all nodes in the phylogeny
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  //calculate intra-sample diversity
  for (plot = 0; plot < S.nsamples; plot++) {
    D[plot] = 0.0;
    Dp[plot] = 0.0;
    for (i = 0; i < S.srec[plot]; i++) {
      for (j = 0; j < S.srec[plot]; j++) {
        if (i != j) {
          D[plot] += S.pabund[plot][i] * S.pabund[plot][j];
          Dp[plot] += (TaxaDist(P, S.id[plot][i], S.id[plot][j],
                                attach) / 2.0) * S.pabund[plot][i]
            * S.pabund[plot][j];
        }
      }
    }
  }

  //calculate inter-sample diversity
  for (k = 0; k < S.nsamples-1; k++) {
    for (l = k+1; l < S.nsamples; l++) {
      H[k][l] = 0.0;
      Hp[k][l] = 0.0;
      for (i = 0; i < S.srec[k]; i++) {
        for (j = 0; j < S.srec[l]; j++) {
          if (S.id[k][i] != S.id[l][j]) {
            H[k][l] += S.pabund[k][i] * S.pabund[l][j];
          }
          Hp[k][l] += (TaxaDist(P, S.id[k][i], S.id[l][j],
                                attach) / 2.0) * S.pabund[k][i]
            * S.pabund[l][j];

        }
      }
      H[k][l] = H[k][l] - (D[k] + D[l])/2.0;
      Hp[k][l] = Hp[k][l] - (Dp[k] + Dp[l])/2.0;
    }
  }

  //calculate total diversity
  Dtotal = 0.0;
  Dtotalp = 0.0;
  for (i = 0; i < S.ntaxa; i++) {
    for (j = 0; j < S.ntaxa; j++) {
      if (i != j) {
        Dtotal += (double) S.psppabund[i] * (double) S.psppabund[j];
        // printf("%d, %d, %f, %f, %f %s\n", i, j, S.psppabund[i], S.psppabund[j], Dtotal, S.taxa[j]);
        Dtotalp += ((double) TaxaDist(P, i, j, attach) / 2.0)
          * (double) S.psppabund[i] * (double) S.psppabund[j];
      }
      sumDist += (TaxaDist(P, i, j, attach) / 2.0);
      countDist++;
    }
  }

  //calculate alpha and beta diversity
  Dalpha = 0.0;
  Dalphap = 0.0;
  for (k = 0; k < S.nsamples; k++) {
    Dalpha += D[k] * ((float)S.irec[k] / (float)S.totabund);
    Dalphap += Dp[k] * ((float)S.irec[k] / (float)S.totabund);
  }
  Dbeta = 0.0;
  Dbetap = 0.0;
  for (k = 0; k < S.nsamples; k++) {
    for (l = 0; l < S.nsamples; l++) {
      if (k != l) {
        Dbeta += H[k][l] * ((float)S.irec[k] / (float)S.totabund)
          * ((float)S.irec[l] / (float)S.totabund);
        Dbetap += Hp[k][l] * ((float)S.irec[k] / (float)S.totabund)
          * ((float)S.irec[l] / (float)S.totabund);
      }
    }
  }

  Dbeta *= 2;
  Dbetap *= 2;

  Fst = Dbeta / Dtotal;
  Fstp = Dbetap / Dtotalp;

  //write results
  if (S.nsamples > 0) {
    for (i = 1; i<S.nsamples; i++) {
      if (S.irec[i] != S.irec[i-1]) {
        printf("#Warning - diversity components are only interpretable when all communities contain equal numbers of individuals. See manual for details.\n");
        break;
      }
    }
  }
  printf("#----- DIVERSITY COMPONENTS\n");
  printf("Component\tStandard\tPhylogenetic\n");
  printf("Alpha\t%f\t%f\n", Dalpha, Dalphap);
  printf("Beta\t%f\t%f\n", Dbeta, Dbetap);
  printf("Total\t%f\t%f\n", Dtotal, Dtotalp);
  printf("Fst\t%f\t%f\n", Fst, Fstp);
  printf("#----- WITHIN-COMMUNITY DIVERSITY\n");
  printf("Plot\tNSpp\tNIndiv\tPropIndiv\tD\tDp\n");
  for (plot = 0; plot < S.nsamples; plot++) {
    printf("%s\t%d\t%d\t%f\t%f\t%f\n", S.pname[plot], S.srec[plot],
           S.irec[plot], (float)S.irec[plot]/(float)S.totabund, D[plot],
           Dp[plot]);
  }

  printf("#----- AMONG-COMMUNITY DIVERSITY (D)\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      if (k == l) {
        printf("\t%f", H[k][l] + (D[k] + D[l])/2.0);
      } else if (k > l) {
        printf("\t%f", H[l][k] + (D[k] + D[l])/2.0);
      } else {
        printf("\t%f", H[k][l] + (D[k] + D[l])/2.0);
      }
    }
    printf("\n");
  }

  printf("#----- AMONG-COMMUNITY DIVERSITY (H)\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      if (k == l) {
        printf("\t%f", H[k][l]);
      } else if (k > l) {
        printf("\t%f", H[l][k]);
      } else {
        printf("\t%f", H[k][l]);
      }
    }
    printf("\n");
  }

  printf("#----- AMONG-COMMUNITY PHYLOGENETIC DIVERSITY (Dp)\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      if (k == l) {
        printf("\t%f", Hp[k][l] + (Dp[k] + Dp[l])/2.0);
      } else if (k > l) {
        printf("\t%f", Hp[l][k] + (Dp[k] + Dp[l])/2.0);
      } else {
        printf("\t%f", Hp[k][l] + (Dp[k] + Dp[l])/2.0);
      }
    }
    printf("\n");
  }

  printf("#----- AMONG-COMMUNITY PHYLOGENETIC DIVERSITY (Hp)\n");
  printf(".");
  for (k = 0; k < S.nsamples; k++)
    printf("\t%s", S.pname[k]);
  printf("\n");
  for (k = 0; k < S.nsamples; k++) {
    printf("%s", S.pname[k]);
    for (l = 0; l < S.nsamples; l++) {
      if (k == l) {
        printf("\t%f", Hp[k][l]);
      } else if (k > l) {
        printf("\t%f", Hp[l][k]);
      } else {
        printf("\t%f", Hp[k][l]);
      }
    }
    printf("\n");
  }

  //TODO free vectors and matrices
  free_dvector(D, 0, S.nsamples-1);
  free_dvector(Dp, 0, S.nsamples-1);
  free_dmatrix(H, 0, S.nsamples-1, 0, S.nsamples-1);
  free_dmatrix(Hp, 0, S.nsamples-1, 0, S.nsamples-1);
  free_ivector(attach, 0, S.ntaxa-1);

}
