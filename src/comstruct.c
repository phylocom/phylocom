// comstruct.c - additional community phylogenetic structure algs

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

// ---------- PD --------------------------------------------------------------
// TODO - modify so usable with randomization methods?
//      - from Ryan Jones (email to C 24may08), calc PD for the virtual
//        community of all taxa in samples

void PD(struct phylo P, struct sample S)
{
  int i, j, x;
  float totalbl = 0.0;
  float *sampletotbl;
  int *nodeVisited;
  int *attach;

  // once only initialize
  sampletotbl = vector(0, S.nsamples-1);
  nodeVisited = ivector(0, P.nnodes-1);
  nodeVisited[0] = 1; // To stop it going further

  attach = ivector(0, S.ntaxa-1);
  AttachSampleToPhylo(S, P, attach);

  // 1. Calc total for the tree:

  // initialize:
  for (j = 1; j < P.nnodes; j++) nodeVisited[j] = 0;

  for (j = 1; j < P.nnodes; j++)
    {
      x = j;
      while (nodeVisited[x] == 0)
        {
          nodeVisited[x] = 1;
          totalbl += P.bl[x];
          x = P.up[x];
        }
    }

  // 2. now for each sample
  for (i=0; i< S.nsamples; i++)
    {
      // initialize
      for (j = 1; j < P.nnodes; j++) nodeVisited[j] = 0;
      sampletotbl[i] = 0.0;

      for (j = 0; j < S.srec[i]; j++)
        {
          x = P.t2n[ attach[ S.id[i][j] ] ];
          while (nodeVisited[x] == 0)
            {
              nodeVisited[x] = 1;
              sampletotbl[i] += P.bl[x];
              x = P.up[x];
            }
        }
    }

  // 3. Output
  printf("sample\tntaxa\tPD\ttreeBL\tpropTreeBL\n");
  for (i=0; i< S.nsamples; i++)
    {
      // "%-7s %7d %7f %7f %7f\n" - for free precision
      //printf("%-7s %7d %7.3f %7.3f %7.3f\n", S.pname[i], S.srec[i], sampletotbl[i],
      //   totalbl, sampletotbl[i]/totalbl);
      printf("%s\t%-d\t%-7.3f\t%-7.3f\t%-7.3f\n", S.pname[i], S.srec[i], sampletotbl[i],\
             totalbl, sampletotbl[i]/totalbl);
    }

}

// ---------- COMDIST --------------------------------------------------------

void ComDist(struct phylo P, struct sample S)
{

  int i, j, p, q, comb;

  float totR;
  float **mdist;
  int *attach; // "make a pointer to an int"

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  AttachSampleToPhylo(S, P, attach);

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  mdist = matrix(0, S.nsamples-1, 0, S.nsamples - 1);

  //intialize output matrix
  for (i= 0; i < S.nsamples; i++)
    {
      for (j= 0; j < S.nsamples; j++)
        {
          mdist[i][j] = 0.0;
        }
    }

  // for each pair of plots:
  for (p=0; p<S.nsamples-1; p++)
    {
      for (q=p+1; q<S.nsamples; q++)
        {

          // for each pair of taxa within each plot
          totR = 0.0;
          comb = 0;

          // all pairs, not pairs within the same set!
          for (i=0; i<S.srec[p]; i++)
            {
              for (j=0; j< S.srec[q]; j++)
                {
                  totR = totR + P.dist[ P.t2n[ attach[ S.id[p][i] ] ] ] [ P.t2n[ attach[ S.id[q][j] ] ] ];
                  comb++;

                }
            }
          mdist[p][q] = totR / (float) comb;
          mdist[q][p] = mdist[p][q];
        }
    }

  // print out result for each plot combo:
  // headers
  printf(".");
  for (i= 0; i < S.nsamples; i++)
    {
      printf("\t%s", S.pname[i]);
    }
  printf("\n");

  // ROWS
  for (i= 0; i < S.nsamples; i++)
    {
      printf("%s", S.pname[i]);

      for (j= 0; j < S.nsamples; j++)
        {
          printf("\t%f", mdist[j][i]);
        }
      printf("\n");
    }

}

// ---------- COMDIST NN ------------------------------------------------------
void ComDistNN(struct phylo P, struct sample S)
{
  int i, j, p, q;
  float minR = 0.0;
  float rel, totminR;

  float **mdist;
  int *attach; // "make a pointer to an int"

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  AttachSampleToPhylo(S, P, attach);

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  mdist = matrix(0, S.nsamples-1, 0, S.nsamples - 1);

  //intialize output matrix
  for (i= 0; i < S.nsamples; i++)
    {
      for (j= 0; j < S.nsamples; j++)
        {
          mdist[i][j] = 0.0;
        }
    }

  // for each pair of plots:
  for (p=0; p<S.nsamples; p++)
    {
      for (q=p+1; q<S.nsamples; q++)
        {

          // for each pair of taxa within each plot
          totminR =  0.0;
          for (i=0; i < S.srec[p]; i++)
            {
              minR = 99999999.0;

              for (j=0; j<S.srec[q]; j++)
                {
                  rel = P.dist[ P.t2n[ attach[ S.id[p][i] ] ] ]\
                    [ P.t2n[ attach[ S.id[q][j] ] ] ];
                  if (minR > rel) {minR = rel;}
                }
              totminR += minR;
            }
          mdist[p][q] = totminR / (float)S.srec[p];
          mdist[q][p] = mdist[p][q];
        }
    }

  // print out result for each plot combo:
  // headers
  printf("SAMPLES");
  for (i= 0; i < S.nsamples; i++)
    {
      printf("\t%s", S.pname[i]);
    }
  printf("\n");

  // ROWS
  for (i= 0; i < S.nsamples; i++)
    {
      printf("%s", S.pname[i]);

      for (j= 0; j < S.nsamples; j++)
        {
          printf("\t%f", mdist[j][i]);
        }
      printf("\n");
    }

}

//-----------------------ICOMDIST-------------------------------

void IComDist(struct phylo P, struct sample S)
{

  int i, j, p, q, comb;
  float totR;
  //float **mdist;
  int *attach; // "make a pointer to an int"
  // float **dmatrix;

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  AttachSampleToPhylo(S, P, attach);

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  // for each pair of plots:
  for (p=0; p<S.nsamples; p++)
    {
      for (q=0; q<S.nsamples; q++)
        {

          // all pairs, not pairs within the same set!
          for (i=0; i<S.srec[p]; i++)
            {

              totR = 0.0;
              comb = 0;

              for (j=0; j < S.srec[q]; j++)
                {
                  totR = totR + P.dist[ P.t2n[ attach[ S.id[p][i] ] ] ][ P.t2n[ attach[S.id[q][j]] ] ];
                  comb++;

                }
              // print results
              printf("AV\t%s\t%s\t%s\t%f\n", S.pname[p], S.taxa[S.id[p][i]], S.pname[q], totR / (float) comb);
            }
        }
    }
  free_matrix(P.dist, 0, P.nnodes-1, 0, P.nnodes-1);
}

//-----------------------ICOMDIST-------------------------------

void IComDistNN(struct phylo P, struct sample S)
{

  int i, j, p, q;
  float minR, R;
  //float **mdist;
  int *attach; // "make a pointer to an int"
  // float **dmatrix;

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  // old: set attach to the pointer passed back by the function -
  //      this is pointless!  It's the same pointer
  // attach = AttachSampleToPhylo(S, P, attach);
  AttachSampleToPhylo(S, P, attach);

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  // for each pair of plots:
  for (p=0; p<S.nsamples; p++)
    {
      for (q=0; q<S.nsamples; q++)
        {

          // all pairs, not pairs within the same set!
          for (i=0; i<S.srec[p]; i++)
            {

              minR = 100000.0;

              for (j=0; j < S.srec[q]; j++)
                {
                  R = P.dist[ P.t2n[ attach[ S.id[p][i] ] ] ][ P.t2n[ attach[S.id[q][j]] ] ];
                  if (R < minR) minR = R;
                }
              // print results
              printf("NT\t%s\t%s\t%s\t%f\n", S.pname[p], S.taxa[S.id[p][i]], S.pname[q], minR);  
            }
        }
    }
  free_matrix(P.dist, 0, P.nnodes-1, 0, P.nnodes-1);
}

//-----------------------VCOMDIST-------------------------------

void VComDist(struct phylo P, struct sample S)
{

  int i, j, p, comb;
  float totR;
  //float **mdist;
  int *attach; // "make a pointer to an int"
  // float **dmatrix;

  means InM = ReadMeans(P, "means");

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  // old: set attach to the pointer passed back by the function -
  //      this is pointless!  It's the same pointer
  // attach = AttachSampleToPhylo(S, P, attach);
  AttachSampleToPhylo(S, P, attach);

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  // for each plot:
  for (p=0; p<S.nsamples; p++)
    {
      // all pairs
      for (i=0; i<S.srec[p]; i++)
        {

          totR = 0.0;
          comb = 0;

          for (j=0; j < S.srec[p]; j++)
            {
              if(i != j)
                {
                  totR = totR + P.dist[ P.t2n[ attach[ S.id[p][i] ] ] ][ P.t2n[ attach[S.id[p][j]] ] ];
                  comb++;
                }
            }
          // print results
          printf("AV\t%s\t%d\t%s\t%f\t%f\n", S.pname[p], S.srec[p], \
                 S.taxa[S.id[p][i]], (totR / (float) comb),  \
                 (-1 * ( (totR / (float) comb) - InM.mpd[S.srec[p]]) / \
                  InM.sdpd[S.srec[p]]));
        }

    }
  free_matrix(P.dist, 0, P.nnodes-1, 0, P.nnodes-1);
}

//-----------------------VCOMDIST-------------------------------

void VComDistNN(struct phylo P, struct sample S)
{

  int i, j, p;
  float minR, R;
  //float **mdist;
  int *attach; // "make a pointer to an int"
  // float **dmatrix;

  means InM = ReadMeans(P, "means");

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  // old: set attach to the pointer passed back by the function -
  //      this is pointless!  It's the same pointer
  // attach = AttachSampleToPhylo(S, P, attach);
  AttachSampleToPhylo(S, P, attach);

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  // for each pair of plots:
  for (p=0; p<S.nsamples; p++)
    {

      // all pairs, not pairs within the same set!
      for (i=0; i<S.srec[p]; i++)
        {

          minR = 100000.0;

          for (j=0; j < S.srec[p]; j++)
            {
              if (i != j)
                {
                  R = P.dist[ P.t2n[ attach[ S.id[p][i] ] ] ][ P.t2n[ attach[S.id[p][j]] ] ];
                  if (R < minR) minR = R;
                }
            }
          // print results
          printf("NN\t%s\t%d\t%s\t%f\t%f\n", S.pname[p], S.srec[p], \
                 S.taxa[S.id[p][i]], minR, \
                 (-1 * ( minR - InM.mnnd[S.srec[p]] ) / \
                  InM.sdnnd[S.srec[p]]));
        }
    }
  free_matrix(P.dist, 0, P.nnodes-1, 0, P.nnodes-1);
}

// ------------------- VMEANS ------------------------------------

void VMeans(struct phylo P)
{

  int i, j, comb, n, r, found;
  float totR, minR, R;
  float mpdTot, mpdTotSqr, mpdTotComb, nndTot, nndTotSqr, nndTotComb;
  float mean, var, NNmean, NNvar;
  int *termvect;
  int *termsamp;
  int *sudah;
  int rnd;

  // allocate a matrix of memory to the pre-existing pointer P.dist
  P.dist = matrix(0, P.nnodes-1, 0, P.nnodes-1);
  DistMatrix(P);

  // make a vector of nodes that correspond to terminal taxa
  termvect = ivector(0, P.termtaxa-1); j = 0;
  for (i=0; i< P.nnodes; i++)
    {
      if (P.noat[i] == 0)
        {
          termvect[j]=i;
          //printf("n%d\tnd%d\n", j, termvect[j]);
          j++;
        }
    }

  // make a sample of nodes that correspond to terminal taxa
  termsamp = ivector(0, P.termtaxa-1);

  // make a list of terminal taxa that have been found
  sudah = ivector(0, P.termtaxa-1);

  // for increasing numbers
  for (n = 2; n < P.termtaxa; n++)
    {
      mpdTot = 0.0;
      mpdTotSqr = 0.0;
      mpdTotComb = 0;
      nndTot = 0.0;
      nndTotSqr = 0.0;
      nndTotComb = 0;

      for (r = 0; r < 1000; r++)
        {
          // make random sample of n taxa from termvect
          found = 0;
          for (i = 0; i < P.termtaxa; i++) sudah[i] = 0;
          while (found < n)
            {
              rnd = (int) (((float) P.termtaxa * (float) random()) / (float) (RAND_MAX+1.0));
              if (sudah[rnd] == 0)
                {
                  termsamp[found] = termvect[rnd];
                  sudah[rnd] = 1;
                  //printf("%d(%d)  ", termsamp[found], rnd);
                  found++;
                }
            }

          // calculate the mean distance from each focal to each other taxon

          for (i=0; i< n; i++)
            {
              totR = 0.0;
              comb = 0;

              for (j=0; j < n; j++)
                {
                  if(i != j)
                    {
                      totR = totR + P.dist[ termsamp[i] ][ termsamp[j] ];
                      comb++;
                    }
                }

              mpdTot += (totR / (float) comb); mpdTotComb++;
              mpdTotSqr += (totR / (float) comb) * (totR / (float) comb);
            }

          // calculate the nn distance to any other taxon in sample
          for (i=0; i<n; i++)
            {

              minR = 100000.0;

              for (j=0; j < n; j++)
                {
                  if (i != j)
                    {
                      R = P.dist[ termsamp[i] ][ termsamp[j] ];
                      if (R < minR) minR = R;
                    }
                }
              nndTot += minR; nndTotComb++;
              nndTotSqr += minR * minR;
            }
        }

      // Calculate stats
      mean = mpdTot / (float) mpdTotComb;
      var   = (mpdTotSqr - ((mpdTot * mpdTot) / (float) mpdTotComb)) \
        / (float) (mpdTotComb - 1);

      // nn
      NNmean = nndTot / (float) nndTotComb;
      NNvar   = (nndTotSqr - ((nndTot * nndTot) / (float) nndTotComb)) \
        / (float) (nndTotComb - 1);

      // There is a slight rounding error that makes levelvar -ve.  Fix:
      // if (var < 0) levelvar *= -1.0;
      // if (NNvar < 0) NNlevelvar *= -1.0;

      // output
      printf("%d\t%f\t%f\t%f\t%f\n", \
             n, (float) mean, (float) sqrt(var), NNmean, sqrt(NNvar));

    }
  free_matrix(P.dist, 0, P.nnodes-1, 0, P.nnodes-1);
}


// ---------- LTT -------------------------------------------------------------
// on the way to being fixed
void Ltt(struct phylo P, struct sample S)
{
  int ltt1 = 0;
  int ltt2 = 0;
  int ltt3 = 0;
  int *deep1;
  int *deep2;
  int *deep3;
  float r2t = 0.0;
  int i, xnode, j;
  float runningd, orunningd;
  int *attach; // "make a pointer to an int"


  deep1 = ivector(0, P.nnodes-1);
  deep2 = ivector(0, P.nnodes-1);
  deep3 = ivector(0, P.nnodes-1);

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"
  AttachSampleToPhylo(S, P, attach);

  // Initialize
  for (i = 0; i < P.nnodes; i++)
    {
      deep1[i] = 0; deep2[i] = 0; deep3[i] = 0;
    }


  xnode = P.t2n[ attach[ S.id[0][0] ] ];
  // Calculate lenghth from root to tips
  while (xnode != 0)
    {
      r2t += P.bl[xnode];
      xnode = P.up[xnode];
    }

  // Ltt for whole phylo
  for (i = 0; i < P.termtaxa; i++)
    {
      xnode = P.t2n[i];
      runningd = 0.0; orunningd = 0.0;

      while (xnode != 0)
        {
          runningd += P.bl[xnode];
          if ((runningd > (r2t * 0.25)) && (orunningd <= (r2t * 0.25)))
            {deep1[xnode] = 1;}
          if ((runningd > (r2t * 0.5)) && (orunningd <= (r2t * 0.5)))
            {deep2[xnode] = 1;}
          if ((runningd > (r2t * 0.75)) && (orunningd <= (r2t * 0.75)))
            {deep3[xnode] = 1;}
          xnode = P.up[xnode];
          orunningd = runningd;
        }
    }

  // Count them up
  for (i = 0; i < P.nnodes; i++)
    {
      if (deep3[i] == 1) {ltt1++;}
      if (deep2[i] == 1) {ltt2++;}
      if (deep1[i] == 1) {ltt3++;}
    }

  // Output
  printf("sample  \ttaxa\tt=0\tt=0.25\tt=0.50\tt=0.75\tt=1.00\n");
  printf("all     \t%d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", P.termtaxa, \
         1.0/(float)P.termtaxa,                     \
         (float) ltt1/(float)P.termtaxa, (float)ltt2/(float)P.termtaxa, \
         (float)ltt3/(float)P.termtaxa, 1.0);

  // repeat for each plot
  for (j = 0; j < S.nsamples; j++)
    {

      // Renitialize
      for (i = 0; i < P.nnodes; i++)
        {
          deep1[i] = 0; deep2[i] = 0; deep3[i] = 0;
        }
      ltt1=0; ltt2 = 0; ltt3=0;

      for (i = 0; i < S.srec[j]; i++)
        {
          xnode = P.t2n[ attach[ S.id[j][i] ] ];
          runningd = 0.0; orunningd = 0.0;

          while (xnode != 0)
            {
              runningd += P.bl[xnode];
              if ((runningd > (r2t * 0.25)) && (orunningd <= (r2t * 0.25)))
                {deep1[xnode] = 1;}
              if ((runningd > (r2t * 0.5)) && (orunningd <= (r2t * 0.5)))
                {deep2[xnode] = 1;}
              if ((runningd > (r2t * 0.75)) && (orunningd <= (r2t * 0.75)))
                {deep3[xnode] = 1;}
              xnode = P.up[xnode];
              orunningd = runningd;
            }
        }

      // Count them up
      for (i = 0; i < P.nnodes; i++)
        {
          if (deep3[i] == 1) {ltt1++;}
          if (deep2[i] == 1) {ltt2++;}
          if (deep1[i] == 1) {ltt3++;}
        }

      // Output
      printf("%-8s\t%d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", S.pname[j], \
             S.srec[j], 1.0/(float)P.termtaxa,(float) ltt1/(float)S.srec[j], \
             (float)ltt2/(float)S.srec[j],  \
             (float) ltt3 / (float) S.srec[j], 1.0);

    }

}

// ---------- LTTR ------------------------------------------------------------

void LttR(struct phylo P, struct sample S)
{
  int run, y, ord1, ord2, ord3, ord4, ord5, ord75;
  int ltt1 = 0;
  int ltt2 = 0;
  int ltt3 = 0;
  int ltt4 = 0;
  int ltt5 = 0;
  int ltt75 = 0;
  float *ltt1_r;
  float *ltt2_r;
  float *ltt3_r;
  float *ltt4_r;
  float *ltt5_r;
  float *ltt75_r;
  int *deep1;
  int *deep2;
  int *deep3;
  int *deep4;
  int *deep5;
  int *deep75;
  float r2t = 0.0;
  int i, xnode, j;
  float runningd, orunningd;
  int *attach; // "make a pointer to an int"
  int   **dummyID;

  attach = ivector(0, S.ntaxa-1); // "allocate memory to a array pointed to
                                  //  by the pointer attach"

  dummyID = imatrix(0, S.nsamples - 1, 0, S.maxrec);

  AttachSampleToPhylo(S, P, attach);

  // Initialize

  ltt1_r = vector(0, RUNS-1);
  ltt2_r = vector(0, RUNS-1);
  ltt3_r = vector(0, RUNS-1);
  ltt4_r = vector(0, RUNS-1);
  ltt5_r = vector(0, RUNS-1);
  ltt75_r = vector(0,RUNS-1);

  deep1 = ivector(0, P.nnodes-1);
  deep2 = ivector(0, P.nnodes-1);
  deep3 = ivector(0, P.nnodes-1);
  deep4 = ivector(0, P.nnodes-1);
  deep5 = ivector(0, P.nnodes-1);
  deep75 = ivector(0, P.nnodes-1);

  for (i = 0; i < P.nnodes; i++)
    {
      deep1[i] = 0; deep2[i] = 0; deep3[i] = 0;
      deep4[i] = 0; deep5[i] = 0; deep75[i] = 0;
    }

  xnode = P.t2n[ attach[ S.id[0][0] ] ];

  // Calculate lenghth from root to tips
  while (xnode != 0)
    {
      r2t += P.bl[xnode];
      xnode = P.up[xnode];
    }

  // Ltt for whole phylo
  for (i = 0; i < P.termtaxa; i++)
    {
      xnode = P.t2n[i];
      runningd = 0.0; orunningd = 0.0;

      while (xnode != 0)
        {
          runningd += P.bl[xnode];
          if ((runningd > (r2t * 0.1)) && (orunningd <= (r2t * 0.1)))
            {deep1[xnode] = 1;}
          if ((runningd > (r2t * 0.2)) && (orunningd <= (r2t * 0.2)))
            {deep2[xnode] = 1;}
          if ((runningd > (r2t * 0.3)) && (orunningd <= (r2t * 0.3)))
            {deep3[xnode] = 1;}
          if ((runningd > (r2t * 0.4)) && (orunningd <= (r2t * 0.4)))
            {deep4[xnode] = 1;}
          if ((runningd > (r2t * 0.5)) && (orunningd <= (r2t * 0.5)))
            {deep5[xnode] = 1;}
          if ((runningd > (r2t * 0.75)) && (orunningd <= (r2t * 0.75)))
            {deep75[xnode] = 1;}
          xnode = P.up[xnode];
          orunningd = runningd;
        }
    }

  // Count them up
  for (i = 0; i < P.nnodes; i++)
    {
      if (deep1[i] == 1) {ltt1++;}
      if (deep2[i] == 1) {ltt2++;}
      if (deep3[i] == 1) {ltt3++;}
      if (deep4[i] == 1) {ltt4++;}
      if (deep5[i] == 1) {ltt5++;}
      if (deep75[i] == 1) {ltt75++;}
    }

  // Output
  printf("sample    taxa\tt=0\tt=0.25\tt=0.50\tt=0.60\tt=0.70\tt=0.80\tt=0.90\tt=1.00\n");
  printf("all       %4d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", P.termtaxa, 1.0/(float)P.termtaxa, (float) ltt75/(float)P.termtaxa, (float)ltt5/(float)P.termtaxa, (float)ltt4/(float)P.termtaxa, (float)ltt3/(float)P.termtaxa, (float)ltt2/(float)P.termtaxa, (float)ltt1/(float)P.termtaxa, 1.0);


  //2010-04-02 major bug fix: The randomizition was scrambling the S.id
  // vector and the 2nd plot was not real!  We use a dummy variable to 
  // store the original data
  for (j = 0; j < S.nsamples; j++) 
    {
      for (i = 0; i < S.srec[j]; i++)
        {
          dummyID[j][i] = S.id[j][i] ;
          // printf("j %d  i %d  id %d\n", j, i, S.id[j][i]);
        }
    }
 

  // repeat for each plot
  for (j = 0; j < S.nsamples; j++)
    {

      // Renitialize
      for (i = 0; i < P.nnodes; i++)
        {
          deep1[i] = 0; deep2[i] = 0; deep3[i] = 0;
          deep4[i] = 0; deep5[i] = 0; deep75[i] = 0;
        }
      ltt1=0; ltt2 = 0; ltt3=0;
      ltt4=0; ltt5 = 0; ltt75=0;

      // reset after each randomization
      AttachSampleToPhylo(S, P, attach);

      // for each species in each sample
      for (i = 0; i < S.srec[j]; i++)
        {
          xnode = P.t2n[ attach[ dummyID[j][i] ] ];
          runningd = 0.0; orunningd = 0.0;
          while (xnode != 0)
            {
              runningd += P.bl[xnode];
              if ((runningd > (r2t * 0.1)) && (orunningd <= (r2t * 0.1)))
                {deep1[xnode] = 1;}
              if ((runningd > (r2t * 0.2)) && (orunningd <= (r2t * 0.2)))
                {deep2[xnode] = 1;}
              if ((runningd > (r2t * 0.3)) && (orunningd <= (r2t * 0.3)))
                {deep3[xnode] = 1;}
              if ((runningd > (r2t * 0.4)) && (orunningd <= (r2t * 0.4)))
                {deep4[xnode] = 1;}
              if ((runningd > (r2t * 0.5)) && (orunningd <= (r2t * 0.5)))
                {deep5[xnode] = 1;}
              if ((runningd > (r2t * 0.75)) && (orunningd <= (r2t * 0.75)))
                {deep75[xnode] = 1;}
              xnode = P.up[xnode];
              orunningd = runningd;
            }
        }

      // Count them up
      for (i = 0; i < P.nnodes; i++)
        {
          if (deep1[i] == 1) {ltt1++;}
          if (deep2[i] == 1) {ltt2++;}
          if (deep3[i] == 1) {ltt3++;}
          if (deep4[i] == 1) {ltt4++;}
          if (deep5[i] == 1) {ltt5++;}
          if (deep75[i] == 1) {ltt75++;}
        }


      // now randomize
      // NOTE!  This is wasteful of randomization runs.  Code written
      // to randomize one plot at a time (RandomizeB(j)), but Steve's null
      // models shuffle all plots.  For speed of fixing, I have not added
      // plot indexes, but maintained code.  Needs fixing.

      for (run = 0; run < RUNS; run++)
        {

          // RandomizeB(j);

          // see comments in combase
          switch (SWAPMETHOD)
            {
              //case 0:
              //PhylogenyAttachShuffle(P, S, attach);
              //break;
            case 1:
              RandomizeSampleTaxaShuffle(S);
              break;
            case 2:
              RandomizeSampleTaxaShuffle(S);
              PhylogenyAttachShuffle(P, S, attach);
              break;
              // case 3:
              // IndependentSwap(S, SWAPS);
              // break;
            default:
              printf("Please use -m command line switch to specify a randomization method.\n");
              printf("See documentation for a list of possible null models.\n");
              exit(EXIT_FAILURE);
              break;
            }

          // check the randomizations
          // for (i = 0; i < S.srec[0]; i++) printf("%d ", S.id[0][i]);
          // printf("\n");


          // Renitialize
          for (i = 0; i < P.nnodes; i++)
            {
              deep1[i] = 0; deep2[i] = 0; deep3[i] = 0;
              deep4[i] = 0; deep5[i] = 0; deep75[i] = 0;
            }
          ltt1_r[run]=0.0; ltt2_r[run]=0.0; ltt3_r[run]=0.0;
          ltt4_r[run]=0.0; ltt5_r[run]=0.0; ltt75_r[run]=0.0;

          for (i = 0; i < S.srec[j]; i++)
            {
              xnode = P.t2n[ attach[ S.id[j][i] ] ];
              // printf("%d ", xnode);
              runningd = 0.0; orunningd = 0.0;
              while (xnode != 0)
                {
                  runningd += P.bl[xnode];
                  if ((runningd > (r2t * 0.1)) && (orunningd <= (r2t * 0.1)))
                    {deep1[xnode] = 1;}
                  if ((runningd > (r2t * 0.2)) && (orunningd <= (r2t * 0.2)))
                    {deep2[xnode] = 1;}
                  if ((runningd > (r2t * 0.3)) && (orunningd <= (r2t * 0.3)))
                    {deep3[xnode] = 1;}
                  if ((runningd > (r2t * 0.4)) && (orunningd <= (r2t * 0.4)))
                    {deep4[xnode] = 1;}
                  if ((runningd > (r2t * 0.5)) && (orunningd <= (r2t * 0.5)))
                    {deep5[xnode] = 1;}
                  if ((runningd > (r2t * 0.75)) && (orunningd <= (r2t * 0.75)))
                    {deep75[xnode] = 1;}
                  xnode = P.up[xnode];
                  orunningd = runningd;
                }
            }

          // Count them up
          for (i = 0; i < P.nnodes; i++)
            {
              if (deep1[i] == 1) {ltt1_r[run]+=1.0;}
              if (deep2[i] == 1) {ltt2_r[run]+=1.0;}
              if (deep3[i] == 1) {ltt3_r[run]+=1.0;}
              if (deep4[i] == 1) {ltt4_r[run]+=1.0;}
              if (deep5[i] == 1) {ltt5_r[run]+=1.0;}
              if (deep75[i] == 1) {ltt75_r[run]+=1.0;}
            }
        }

      // now the Sorts:

      Sort(ltt1_r, RUNS);
      ord1 = 0;
      for (y = 0; y < RUNS; y++)
        {
          if ((float) ltt1 < ltt1_r[y]) { ord1 = y; break; }
        }
      if (y == RUNS) ord1 = RUNS;  // not stopped

      Sort(ltt2_r, RUNS);
      ord2 = 0;
      for (y = 0; y < RUNS; y++)
        {
          if ((float)ltt2 < ltt2_r[y]) { ord2 = y; break; }
        }
      if (y == RUNS) ord2 = RUNS;  // not stopped

      Sort(ltt3_r, RUNS);
      ord3 = 0;
      for (y = 0; y < RUNS; y++)
        {
          if ((float)ltt3 < ltt3_r[y]) { ord3 = y; break; }
        }
      if (y == RUNS) ord3 = RUNS;  // not stopped

      Sort(ltt4_r, RUNS);
      ord4 = 0;
      for (y = 0; y < RUNS; y++)
        {
          if ((float)ltt4 < ltt4_r[y]) { ord4 = y; break; }
        }
      if (y == RUNS) ord4 = RUNS; // not stopped

      Sort(ltt5_r, RUNS);
      ord5 = 0;
      for (y = 0; y < RUNS; y++)
        {
          if ((float)ltt5 < ltt5_r[y]) { ord5 = y; break; }
        }
      if (y == RUNS) ord5 = RUNS; // not stopped

      Sort(ltt75_r, RUNS);
      ord75 = 0;
      for (y = 0; y < RUNS; y++)
        {
          if ((float)ltt75 < ltt75_r[y]) { ord75 = y; break; }
        }
      if (y == RUNS) ord75 = RUNS; // not stopped


      // Output
      printf("%-8s  %4d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", S.pname[j], S.srec[j], 1.0/(float)S.srec[j], (float) ltt75/(float)S.srec[j], (float)ltt5/(float)S.srec[j], (float)ltt4/(float)S.srec[j], (float)ltt3/(float)S.srec[j], (float)ltt2/(float)S.srec[j], (float)ltt1/(float)S.srec[j], 1.0);
      printf("                      \t%d\t%d\t%d\t%d\t%d\t%d\t\n", ord75, ord5, ord4, ord3, ord2, ord1);

    }

}

void Sort(float Array[], int Nia)
{
  int i;
  float swap;
  int aSwap = 1;

  while (aSwap == 1)
    {
      aSwap = 0;

      for(i = 0; i < Nia -1; i++) /* note up to one less than final Total... */
        {
          if (Array[i] > Array[i+1])
            {
              swap = Array[i+1];
              Array[i+1] = Array[i];
              Array[i] = swap;
              aSwap = 1;
            }
        }
    }
}

//// ---------- SLIDE -----------------------------------------------------------
//
//void Slide()
//{
//
//  int y, run, ordT, plot, w, x;
//
//  float VR_pw[MAXPLOTS][50];
//  float VR_rw[RUNS][50];
//  float VR_r[RUNS];
//
//  int maxwidth;
//
//  // printf(">>> %f\n", SlidingN(0, 5));
//
//  printf("plt\tntx\twidth\tVWR\tv_gt\truns\n");
//
//  // Get results for real data:
//
//  for (plot = 0; plot < Plots; plot++)
//    {
//      if (Plottot[plot]<50) {maxwidth = Plottot[plot];} else {maxwidth = 50;};
//      for (w = 1; w < maxwidth; w++)
//  {
//    VR_pw[plot][w] = SlidingN(plot, w);
//    printf("%d\t%d\t%d\t%f\t1000\n", plot, Plottot[plot],w, VR_pw[plot][w]);
//  }
//    }
//
//  //exit(EXIT_SUCCESS);
//
//  // Now the randomizations:
//  for (plot = 0; plot < Plots; plot++)
//    {
//      if (Plottot[plot]<50) {maxwidth = Plottot[plot];} else {maxwidth = 50;};
//      for (run = 0; run < RUNS; run++)
//        {
//          RandomizeB(plot);
//
//
//    for (w = 1; w < maxwidth; w++)
//      {
//        VR_rw[run][w] = SlidingN(plot, w);
//      }
//  }
//
//
//       NV order
//      //pass the data
//      for (w = 1; w < maxwidth; w++)
//  {
//    for(x = 0; x < RUNS; x++) {VR_r[x] = VR_rw[x][w];};
//    Sort(VR_r, RUNS);
//    ordT = 0;
//    for (y = 0; y < RUNS; y++)
//      {
//        if (VR_pw[plot][w] < VR_r[y])
//      {
//        ordT = y;
//        break;
//      }
//      }
//
//    if (y == RUNS) ordT = RUNS;  not stopped
//
//    printf("%3d\t%3d\t%3d\t%5.2f\t%4d\t%4d\n", \//
//             plot+1, Plottot[plot], w, VR_pw[plot][w], RUNS - ordT, RUNS);
//  }
//    }
//}

//// ---------- SLIDING N -------------------------------------------------------
//
//float SlidingN(int plot, int width)
//{
//  int j, k, l, memb = 0;
//  float rel[MAXINPLOT];
//  float slidecon = 0.0;
//
//  /* for every x */
//  for (j = 0; j < Plottot[plot]; j++)
//    {
//      memb = 0;
//      /* there is a pair with y */
//      for (k = 0; k < Plottot[plot]; k++)
//        {
//          if (k != j)
//            {
//              rel[memb] = Relatedness(Id[plot][j], Id[plot][k]);
//              memb++;
//            }
//        }
//
//      Sort(rel, memb);
//      for (l = 0; l < width; l++)
//  {
//    // printf("%d %d %f\n", j, l, rel[l]);
//    slidecon += rel[l];
//  }
//    }
//  return ((slidecon / (float) width) / (float) Plottot[plot]);
//}

//// ---------- NODE SIG --------------------------------------------------------

void NodeSig(phylo P, sample S, int outmethod, int abundWeighted) {
  // Currently need to use taxon name of interior as a marker,
  // because Mesquite will not show all labels, but if it will soon,
  // best to use notes, and not muck around with taxon name

  //TODO modifying to use abundances, need longs instead of ints for counters

  int plot, node, taxon, i, ordHI, ordLO, run;
  int tipsReal_n[P.nnodes];
  float test_r[RUNS];
  char mark[2];
  char tmp[15];
  int **tips_rn;
  // was: int tips_rn[RUNS][MaxNode+1];
  phylo Out[S.nsamples];
  int *attach;

  attach = ivector(0, S.ntaxa-1);

  tips_rn = imatrix(0, RUNS-1, 0, P.nnodes-1);

  for (i = 0; i < S.nsamples; i++) {
    Out[i] = P; // all the pointers in Out are the same as those in Intree
                // careful not to change any preexisting arrays in Out, or they will
                // also be changed in Intree!

    if (TreeView == 0)
      Out[i].arenotes = 1;
    if (TreeView == 1)
      Out[i].arenotes = 0;

    // dimension a new array for Out names:
    if (TreeView == 0)
      Out[i].notes = cmatrix(0, P.nnodes-1, 0, MAXNOTELENGTH+10);
    if (TreeView == 1)
      Out[i].taxon = cmatrix(0, P.nnodes-1, 0, MAXTAXONLENGTH+10);

  }

  // for each plot

  if (outmethod == 1)
    printf("plot\tnode\tnode_name           \tntaxa\tmedian\trank\tsig\n");
  for (plot = 0; plot < S.nsamples; plot++) {
    if (S.srec[plot] > 2) {

      for (node = 0; node < P.nnodes; node++) {
        tipsReal_n[node] = 0;
        for (run = 0; run < RUNS; run++)
          tips_rn[run][node] = 0;
      }

      // need to reset it
      AttachSampleToPhylo(S, P, attach);

      // follow up from tips, adding 1 to each node passed through
      for (taxon = 0; taxon < S.srec[plot]; taxon++) {
        i = P.t2n[ attach[ S.id[plot][taxon] ] ];
        while (i != -1) {
          if (abundWeighted)
            tipsReal_n[i] += S.abund[plot][taxon];
          else
            tipsReal_n[i]++;
          i = P.up[i];
        }
      }
      // ooo this is slow, putting rnd inside the plot loop!
      for (run = 0; run < RUNS; run++) {
        // now randomize the plot
        // RandomizeB(plot);
        PhylogenyAttachShuffle(P, S, attach);

        for (taxon = 0; taxon < S.srec[plot]; taxon++) {
          i = P.t2n[ attach[ S.id[plot][taxon] ] ];
          while (i != -1) {
            if (abundWeighted)
              tips_rn[run][i] += S.abund[plot][taxon];
            else
              tips_rn[run][i]++;
            i = P.up[i];
          }
        }
      }

      // now unpeel the nodes
      for (node = 0; node < P.nnodes; node++) {
        if (TreeView == 0)
          strcpy(Out[plot].notes[node], "");
        if (TreeView == 1)
          strcpy(Out[plot].taxon[node], "");

        // interior nodes only
        if (P.noat[node] != 0) {
          ordHI = 0;
          ordLO = 0;
          for (run = 0; run < RUNS; run++) {
            if (tips_rn[run][node] < tipsReal_n[node])
              ordHI++;
            if (tips_rn[run][node] > tipsReal_n[node])
              ordLO++;
          }

          strcpy(mark, " ");

          if (ordLO >= (int) ((float) RUNS * 0.975)) {
            strcpy(mark, "-");
            if (TreeView == 0)
              strcat(Out[plot].notes[node], "SIGLESS");
            if (TreeView == 1) {
              sprintf(tmp, "LESS_%d_", node);
              strcpy(Out[plot].taxon[node], tmp);
            }
            if (outmethod == 1)
              printf("%d\t%d\t%-10s\t%d\t%d\t%d\t%d\t%s\n", plot
                     +1, node, P.taxon[node], tipsReal_n[node],
                     (int) test_r[(int) ((float) RUNS * 0.5)], ordHI,
                     ordLO, mark);
          } else if (ordHI >= (int) ((float) RUNS * 0.975)) {
            strcpy(mark, "+");
            if (TreeView == 0)
              strcat(Out[plot].notes[node], "SIGMORE");

            if (TreeView == 1) {
              sprintf(tmp, "MORE_%d_", node);
              strcpy(Out[plot].taxon[node], tmp);
            }
            if (outmethod == 1)
              printf("%d\t%d\t%-20s\t%d\t%d\t%d\t%d\t%s\n", plot
                     +1, node, P.taxon[node], tipsReal_n[node],
                     (int) test_r[(int) ((float) RUNS * 0.5)], ordHI,
                     ordLO, mark);
          }
          // printf("%s\t%s\n", Out.taxon[node], Intree.taxon[node]);
          else if (outmethod == 1)
            printf("%d\t%d\t%-20s\t%d\t%d\t%d\t%d\t%s\n", plot+1,
                   node, P.taxon[node], tipsReal_n[node],
                   (int) test_r[(int) ((float) RUNS * 0.5)], ordHI,
                   ordLO, mark);
        }
        if ((TreeView == 1) && (strcmp(P.taxon[node], ".") != 0)) {
          strcat(Out[plot].taxon[node], P.taxon[node]);
        }
      }

    }
    // Name tree
    strcpy(Out[plot].phyname, "NodeSig_");
    strcat(Out[plot].phyname, S.pname[plot]);
  }
  if (outmethod ==0)
    WriteNexus(Out, S.nsamples, ReadSample(SampleFile), 1,
               ReadTraits(TraitFile), 1);

}

