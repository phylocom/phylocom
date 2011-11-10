#include <stdio.h>  
#include <stdlib.h>  
#include <math.h>  
#include <string.h>  
#include "phylocom.h"  
#include "nrutil.h" 
#include "stats.h" 

// ------ New AOT ----------------------------------------------- //
// New aot function for version 4, all broken into modules, etc.

void AOT(traits T, phylo P, int Output, int xVar)
{
  if (Debug) printf("enter aot\n");

  int i, t, run;
  nodes N, N0;
  int *PTattach;
 
  PTattach = ivector(0,P.nnodes);
  AttachPhyloToTraits(P,T,PTattach);
   
  P.upo = ivector(0,P.nnodes-1); // up pass order
  P.tbl = vector(0,P.nnodes-1); // transformed branch lengths
  P.ntip = ivector(0,P.nnodes-1); // number of tips above node
  P.nint = ivector(0,P.nnodes-1); // number of internals above node
  P.tiplist = imatrix(0,P.nnodes-1,0,P.nnodes-1); // list of tips above node
  P.intlist = imatrix(0,P.nnodes-1,0,P.nnodes-1); // list of internals ""
  
  N0.ordTrt = ivector(0,T.ntraits);

  // tCh, index3: 0 = tipmean; 1 = tipsd; 2 = tipcorr; 
  N.tCh = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);
  N0.tCh = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);
  N0.tChLSig = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);
  N0.tChHSig = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);

  // nCh, index 3: 0 = node mean, wtd by BL; 1 = nodeSD; 2 = piccorr;
  N.nCh = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);
  N0.nCh = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);
  N0.nChLSig = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);
  N0.nChHSig = f3tensor(0,P.nnodes-1,0,T.ntraits-1,0,2);

  // iCon: contrasts for each trait (based on xVar as independent)
  N0.iCon = matrix(0,P.nnodes-1,0,T.ntraits-1);
 
  // cSt: contrast standard deviation
  //      (same for all traits, NA if x trait is binary) 
  N0.cSt = vector(0,P.nnodes-1);
 
  if (Debug) printf("aot stop 0\n");
  
  // rndArr: vector to hold random tip sorts
  N.rndArr = imatrix(0,P.termtaxa-1,0,T.ntraits-1);
  N0.rndArr = imatrix(0,P.termtaxa-1,0,T.ntraits-1);
  
  for (i = 0; i < P.termtaxa; i++) 
    for (t = 0; t < T.ntraits; t++) 
      {
        N0.rndArr[i][t] = i;
        //printf("%d\t%d\n", i, N.rndArr[i][t]);
      }
    
  if (Debug) printf("aot stop 1\n");
  
  // Set trait order (-x run switch)
  N0.ordTrt[0] = xVar - 1;
  t = 0;
  for (i = 1; i < T.ntraits; i++)
    {
      if (t == xVar-1) t++;
      N0.ordTrt[i] = t;
      t++;
    }

  MakeUpPassOrder(P);
  AssignNodeLists(P);
  if (Debug) printf("aot stop 2\n");
  
  // Single character analyses for observed data
  TipStats(P,T,N0,PTattach);
  NodeCharF(P,T,N0,PTattach);
  //AOV(P,T,N); // Node variance analysis
    
  if (Debug) printf("aot stop 3\n");    
  for (run = 1; run <= RUNS; run++)
    {
      RandArray(run, T.ntaxa, T.ntraits, N); 
      TipStats(P,T,N,PTattach);
      NodeCharF(P,T,N,PTattach);
      //AOV(P,T,N); // Node variance analysis
      SigCount(P.nnodes, T.ntraits, N, N0);
    }

  // Pairwise character analyses
  if (T.type[xVar] == 0) binPIC(P,T,N0);
  if (T.type[xVar] == 3) PIC(P,N0);

  aot_outfile(P, T, N0);
}

float *summaryStats(int N,float *X)
{
  int i;
  float sum, sumsq;
  float *stats;

  stats = vector(0,1);
  sum = 0;
  sumsq = 0;
  for (i = 0; i < N; i++)
    {
      sum += X[i];
      sumsq += X[i]*X[i];
    }
  stats[0] = sum/N;
  stats[1] = (sumsq-(sum*sum)/N)/(N-1);
  if (stats[1] < 0) stats[1] = 0;
  stats[1] = sqrt(stats[1]);

  return(stats);
}

float correlation(int N,float *X,float *Y)
{
  int i;
  float R;
  
  R = 0;
  for (i = 0; i < N; i++)
    {
      R += X[i]*Y[i];
    }
  return(R);
}

float *TraitsAtNode(int t, int node, phylo P, traits T, nodes N, int *PTattach)
{
  int i,tmp;
  float *X;
   
  X = vector(0, P.termtaxa-1);
  for (i = 0; i < P.ntip[node]; i++)
    {
      tmp = P.tiplist[node][i];
      X[i] = T.tr[N.rndArr[PTattach[tmp]][t]][t];
    }
  return(X);
  free_vector(X,0,P.termtaxa-1);
}

void NodeCharF(phylo P, traits T, nodes N, int *PTattach)
{
  // Calculation of internal node values based on
  // Felsenstein's algorithm
  int n,j,t,node,upn;
  float ***tmp;
  //float tbl, dsum, dsumsq, dN, tval, nval;
  float tbl, tval;
  tval = 0.0;
  
  tmp = f3tensor(0,P.nnodes-1,0,T.ntraits,0,1);
  
  //printf("Enter InternalTraitValues\n");
  for (n = 0; n < P.nnodes; n++)
    for (t = 0; t < T.ntraits; t++)
      for (j = 0; j <= 1; j++)
        {
          tmp[n][t][j] = 0;
        }

  for (n = 0; n < P.nnodes; n++)
    {
      node = P.upo[n];
      upn = P.up[node];
      if (NoBL == 0) 
        {
          tbl = P.bl[node]; 
        }
      else 
        {
          tbl = 1;
        }
      
      for (t = 0; t < T.ntraits; t++)
        //<<<<<<< aot.c
        //=======
        {
          if (P.noat[node] == 0) //Terminal
            // Terminals
            {
              // NEED TO PUT rndArr back in
              N.nCh[node][t][0] = T.tr[PTattach[node]][t];
              tmp[upn][t][0] += T.tr[PTattach[node]][t] * tbl;
              tmp[upn][t][1] += 1 / tbl;
              P.tbl[node] = P.bl[node];
            }
          //else // NoBL = 1, use equal BL
          //>>>>>>> 3.5
          //{
          //<<<<<<< aot.c
          //tval = T.tr[N.rndArr[T.p2t[node]][t]][t];
          //if (P.noat[node] == 0) //Terminal
          // Terminals
          //=======
          //tmp[upn][t][0] += T.tr[PTattach[node]][t];
          //P.tbl[node] = 1;
          //tmp[upn][t][1]++;
          //}
          //}
          else
            // Internals
            {
              N.nCh[node][t][0] = tmp[node][t][0] /
                tmp[node][t][1];
              //if (t == 0)
              //{
              //if (NoBL == 0) 
              //>>>>>>> 3.5
              //{
              N.nCh[node][t][0] = tval;
              tmp[upn][t][0] += tval * (1/tbl);
              tmp[upn][t][1] += 1 / tbl;
              P.tbl[node] = tbl;
              //}
              //else
              // Internals
              //{
              //N.nCh[node][t][0] = tmp[node][t][0] /
              //        tmp[node][t][1];
              //      if (t == 0)
              //        {
              //          P.tbl[node] = tbl +
              //            (P.noat[node]/tmp[node][t][1]);
              //        }
              // Calculate node standard deviation
              //    d1 = 0;
              //    while (P.up[P.upo[d1]] != node) d1++; //finds first down node
            
              //d1 = P.ldown[node];
              //nval = N.nCh[d1][t][0];
              //dsum = nval;
              //dsumsq = nval * nval;
              //dN = 1;
              //while (P.rsister[d1] != -99)
              //{
              //d1 = P.rsister[d1];
              //nval = N.nCh[d1][t][0];
              //dsum += nval;
              //dsumsq += nval * nval;
              //dN++;
              //}
              //if (dN != P.noat[node]) printf("daughter number mismatch");
              //    N.nCh[node][t][1] = (dsumsq - (dsum*dsum)/dN) /
              //(dN - 1);
              // end node sd calc
            
              if (FALSE) printf("%d\t%d\t%f\t%f\n",
                                node, t, N.nCh[node][t][0], P.tbl[node]);
              if (node != 0)
                {
                  tmp[upn][t][0] += N.nCh[node][t][0] /
                    tbl;
                  tmp[upn][t][1] += 1 / tbl;
                }
            }
        }
    }
}

void TipStats(phylo P, traits T, nodes N, int *PTattach)
{
  int i,j,k,x,y,node;
  float RXY; 
  float *stats;
  float *X;
  float *Y;

  stats = vector(0,1);
  X = vector(0, P.termtaxa-1);
  Y = vector(0, P.termtaxa-1);
  for (node = 0; node < P.nnodes; node++)
    {
      for (j = 0; j < T.ntraits; j++)
        {
          x = N.ordTrt[j];
          X = TraitsAtNode(j,node,P,T,N,PTattach); 
          stats = summaryStats(P.ntip[node], X);
          for (i = 0; i <= 1; i++) 
            {
              //printf("%d\t%f\n", i, stats[i]);
              N.tCh[node][j][i] = stats[i];
            }

          for (k = j+1; k < T.ntraits; k++)
            {
              y = N.ordTrt[k];
              if (P.noat[node] != 0)
                {
                  //TraitsAtNode(trait#, node#, trait struct, phylo struct)
                  Y = TraitsAtNode(k,node,P,T,N,PTattach);
                  RXY = correlation(P.ntip[node],X,Y);
                  //printf("%d\t%d\t%d\t%f\n", j,k,node,RXY);
                }
            }
        }
    }

  if (FALSE)
    {
      for (node = 0; node < P.nnodes; node++)
        for (x = 0; x < T.ntraits; x++)
          {
            printf("%d\t%d\t%f\t%f\n", node, x, N.tCh[node][x][0],N.tCh[node][x][1]);
          }
    } 
  // To HERE

}

void SigCount(int nnodes, int ntraits, nodes N, nodes N0)
{
  int i, node, t;
  for (node = 0; node < nnodes; node++)
    for (t = 0; t < ntraits; t++)
      for (i = 0; i < 2; i++)
        {
          //printf("%f\t%f\n", N0.tCh[node][t][i], N.tCh[node][t][i]);
          if (N0.tCh[node][t][i] <= N.tCh[node][t][i]) \
            N0.tChLSig[node][t][i]++;
          if (N0.tCh[node][t][i] >= N.tCh[node][t][i]) \
            N0.tChHSig[node][t][i]++;
        }
}

void PIC(phylo P, nodes N)
{
  int node, upn;


  for (node = 0; node < P.nnodes; node++)
    {
      upn = 0;

    }
}

void binPIC(phylo P, traits T, nodes N)
{
  int t;

  t = 0;
  if (t == 0)
    {
    }
}

void aot_outfile(phylo P, traits T, nodes N0)
{
  int i, node, t, n;
  for (n = 0; n < P.nnodes; n++)
    {
        
      node = P.nnodes - n - 1;
      if (P.noat[node] != 0)
        {
          for (t = 0; t < T.ntraits; t++)
            {
              printf("%d\t%d\t%d\t", node, t, P.ntip[node]);
              for (i = 0; i < 2; i++)
                {
                  printf("%#1.3f\t%#4.f\t%#4.f\t", N0.tCh[node][t][i], \
                         N0.tChLSig[node][t][i], N0.tChHSig[node][t][i]);
                }
              printf("\n");
            }
        }
    }
}

void RandArray(int run, int nTaxa, int nTraits, nodes N)
{
  // Produces a random list from 0 to n-1 inclusive for each trait
  // IF TRAITS ARE KEPT TOGETHER, THEN TEST IS FOR EFFECTS OF PHYLOGENY
  int data[nTaxa];
  int i, j, loop, size, pick, t;
    
  for (t = 0; t < nTraits; t++)
    {
      for (i = 0; i < nTaxa; i++)
        {
          data[i] = i;
        }
      for (loop = 0; loop < nTaxa; loop++)
        {
          size = nTaxa - loop;
          // pick = ((int) (((float) size * random()) / (RAND_MAX+1.0)));
          pick = (int)((((double) random ())/((double) RAND_MAX))*(double) size);
          if (pick == size)
            {
              pick = pick - 1;
            }
          if (pick > size)
            {
              printf ("RndArr out of bounds, run = %d, loop = %d", run,loop);
            }
          N.rndArr[loop][t] = data[pick];
          for (j = pick; j < size - 1; j++)
            {
              data[j] = data[j + 1];
            }
          if (Debug) printf("%d\t%d\n", t, N.rndArr[loop][t]);
          //printf("%4d %4d %7.3f\n",loop,N.rndArr[loop][t],InC.tr[N.rndArr[loop][t]][t]);
        }
    }
}

