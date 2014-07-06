// prune.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

struct phylo Prune(struct phylo A, int *keep) // array keep is nnodes long
{
  //  phylo A;
  phylo B;
  int *yes;
  int *newup, *ldown, *rsister, *fixednode;
  float *newbl;
  int x, y, mark, newnnodes, newn, d, move, newx, newy;
  // int keep[8] = { 0, 1, 1, 0, 0, 0, 1, 1 };


  // test: set up A -------
  /*    A.nnodes = 8; */
  /*    A.up = ivector(0, A.nnodes-1); */
  /*    A.taxon = cmatrix(0, A.nnodes-1, 0, 10); */
  /*    A.bl = vector(0, A.nnodes-1); */

  /*    A.up[0] = -1; A.up[1] = 0; A.up[2] = 0; A.up[3] = 0;  */
  /*    A.up[4] = 3; A.up[5] = 3; A.up[6] = 5; A.up[7] = 5;  */
  /*    strcpy(A.taxon[0], "root"); strcpy(A.taxon[1], "sp1"); */
  /*    strcpy(A.taxon[2], "sp2"); strcpy(A.taxon[3], "."); */
  /*    strcpy(A.taxon[4], "sp4"); strcpy(A.taxon[5], "."); */
  /*    strcpy(A.taxon[6], "sp6"); strcpy(A.taxon[7], "sp7"); */
  /*    A.bl[0] = 1.0; A.bl[1] = 3.0; A.bl[2] = 3.0; A.bl[3] = 1.0;  */
  /*    A.bl[4] = 2.0; A.bl[5] = 1.0; A.bl[6] = 1.0; A.bl[7] = 1.0;  */
  // -----------------------

  // check:

  y = 0; mark = 0;
  for (x = 0; x < A.nnodes; x++) 
    {
      y += keep[x];
      if (y == 1) mark = x;
    }
  if (y == 0) 
    {
      fprintf(stderr, "Error: No taxa to keep\n");
      exit(8);
    }
  if (y == 1) 
    {
      fprintf(stderr, "Error: one taxon only\n");
      exit(8);
      //      printf("(%s);\n", A.taxon[mark]);
      // return 1;
    }

  // Dimension:
  yes = ivector(0, A.nnodes-1);
  newup = ivector(0, A.nnodes-1);
  ldown = ivector(0, A.nnodes-1);
  rsister = ivector(0, A.nnodes-1);
  newbl =  vector(0, A.nnodes-1);
  fixednode = ivector(0, A.nnodes-1);

  for (x = 0; x < A.nnodes; x++)
    {
      yes[x] = 0;
      newup[x] = -99;
      ldown[x] = -99;
      rsister[x] = -99;
      newbl[x] = 0.0;
      fixednode[x] = 0;
    }


  // run up tree, counting how many live nodes there are
  for (x = 0; x < A.nnodes; x++)
    {
      //printf("%d\t%d\t%d\t%s\n",keep[x], x, A.up[x], A.taxon[x]);
      if (keep[x] == 1)
        {
          yes[x]++;
          y = x;
          while (y != 0)
            {
              y = A.up[y];
              yes[y]++;
            }
        }
    }

  newnnodes = 0;
  // again, making a new structure
  for (x = 0; x < A.nnodes; x++)
    {

      // working root wards from `keepers'
      if (keep[x] == 1)
        {
          newnnodes++; // once for each tip
          y = x;
          // while not yet at the root and not retracing steps
          while ( (y != 0) && (newup[y] == -99) )
            {
              // move up the placeholder 
              newup[y] = A.up[y];
              newbl[y] = A.bl[y];
              // does the new node have 1 descendent?
              while ((yes[newup[y]] == yes[y]) && (newup[y] != 0) )
                {
                  // if so, move up again
                  newbl[y] += A.bl[newup[y]];
                  newup[y] = A.up[newup[y]];

                }
              // now create a dwonwards looking structure
          
              // is this the first daughter in new structure?
              if (ldown[newup[y]] == -99) 
                {
                  ldown[newup[y]] = y;
                  newnnodes++; // once for each new interior node 
                }
              // if not, find the dangling sister
              else if (y != ldown[newup[y]])
                {
                  mark = ldown[newup[y]];
                  //printf("  y%d newup%d ldown%d\n", y, newup[y], mark);
                  while (rsister[mark] != -99)
                    {
                      //printf("                     move%d\n", mark);
                      mark = rsister[mark];
                    }
                  if (mark != y) 
                    {
                      rsister[mark] = y;
                      //printf("                     mark%d to rsis%d\n", mark, y);
                    }
                  // need this check because if a node is passed through twice
                  // it will be added twice
                }
              y = newup[y];
            }
        }
    }
  //    printf("node  __up  -yes  newu  __bl  newb  ldau  rsis\n");
  //for (x = 0; x < A.nnodes; x++)
  //  {
  //    printf("%4d  %4d  %4d  %4d  %4.1f  %4.1f  %4d  %4d  %s\n", x, A.up[x], yes[x], newup[x], A.bl[x], newbl[x], ldown[x], rsister[x], A.taxon[x]);
  //  }
 
  // new tree
  B.nnodes = newnnodes;
  B.up = ivector(0, newnnodes-1);
  B.noat = ivector(0, newnnodes-1);
  B.taxon = cmatrix(0, newnnodes-1, 0, MAXTAXONLENGTH);
  B.bl = vector(0, newnnodes-1);
  B.rsister = ivector(0, newnnodes-1);
  B.ldown = ivector(0, newnnodes-1);
  B.depth = ivector(0, newnnodes-1);
  B.arenotes = 0;

  for (x = 0; x < B.nnodes; x++) B.noat[x] = 0;

  // now, work down:

  newn = 0;
  B.bl[0] = A.bl[0];
  B.up[0] = -1;
  strcpy(B.taxon[0], A.taxon[0]);
  B.depth[0] = 0;
  d = 0;
  x = y = 0;
  newx = newy = 0;
  move = 1;
  
  while (move != -2) // the warning that all possibilities are finished
    {
      // down
      if (move == 1)
        {
          newn++;
          d++;
          y = ldown[x];
          newy = newn;
          B.ldown[newx] = newy;
          B.up[newy] = newx;
          B.depth[newy] = d;
          B.bl[newy] = newbl[y];
          B.noat[newx]++;
          strcpy(B.taxon[newy], A.taxon[y]);

          if (ldown[y] != -99) move = 1;
          else move = 0;
          //printf("x=%d  y=%d  d=%d  newn=%d  nat%d++ moving down\n", x, y, d, newy, newx);
          x = y;
          newx = newy;
        }

      else if (move == 0)
        {
          newn++;
          newy = newn;
          y = rsister[x];
          B.rsister[newx] = newy;
          B.up[newy] = B.up[newx];
          B.depth[newy] = d;
          B.bl[newy] = newbl[y];
          B.noat[B.up[newy]]++;
          strcpy(B.taxon[newy], A.taxon[y]);

          if (ldown[y] != -99) move = 1;
          else if (rsister[y] != -99) move = 0;
          else move= -1;
          //printf("x=%d  y=%d  d=%d  newn=%d  nat%d++ moving right\n", x, y, d, newn, B.up[newy]);
          x = y;
          newx = newy;
        }
      else if (move == -1)
        {
          y = newup[x];
          newy = B.up[newx];
          d--;
          if (rsister[y] != -99) move = 0;
          else if (y != 0) move = -1;
          else move = -2;
          //printf("x=%d  y=%d  d=%d  newn=%d         moving up\n", x, y, d, newn);
          x = y;
          newx = newy;
        }
    }

  //printf("node  __up  noat  __bl  dept  taxon\n");
  //for (x = 0; x < B.nnodes; x++)
  // {
  //    printf("%4d  %4d  %4d  %4.1f  %4d  %s\n", x, B.up[x], B.noat[x], B.bl[x], B.depth[x], B.taxon[x]);
  //  }

  //Fy2newRec(B);  // works

  free_ivector(yes, 0, A.nnodes-1);
  free_ivector(newup,0, A.nnodes-1);
  free_ivector(ldown, 0, A.nnodes-1);
  free_ivector(rsister, 0, A.nnodes-1);
  free_vector(newbl, 0, A.nnodes-1);

  return B;

}

void RandPrune(struct phylo P, int ntx, int ntre)
{
  // make a dummy sample
  struct sample S;
  int i, j;
  int *attach;
  int *keep;

  if (ntx >= P.termtaxa) {
    fprintf(stderr, "Too many taxa for phylo\n"); exit(EXIT_FAILURE);
  }
  
  keep = ivector(0,P.nnodes-1);

  // define some aspects of the sample - not sure if they're all needed, but 
  // better safe 
  S.nsamples = 1;
  S.nrec = ntx; 
  S.maxrec = ntx;
  S.pname = cmatrix(0,0,0,10);
  strcpy(S.pname[0], "shuffle");
  S.srec = ivector(0,0);
  S.srec[0] = ntx;
  S.id = imatrix(0,0,0,ntx-1);
  S.abund = imatrix(0,0,0,ntx-1);
  S.taxa = cmatrix(0,ntx-1,0,MAXTAXONLENGTH);
  for (i = 0; i < ntx; i++) {
    S.id[0][i] = i; S.abund[0][i] = 1;
    strcpy(S.taxa[i], P.taxalist[i]);
  }
  S.ntaxa = ntx;

  attach = ivector(0, ntx-1);
  AttachSampleToPhylo(S, P, attach);
  
  //for (i = 0; i < S.srec[0]; i++) \//
  //    printf("  %-4s", P.taxalist[attach[S.id[0][i]]]);
  //printf("\n");
  
  //for (i = 0; i < P.nnodes; i++) keep[i] = 0;
  //for (i = 0; i < S.srec[0]; i++) keep[P.t2n[attach[S.id[0][i]]]] = 1;
  //for (i = 0; i < P.nnodes; i++) printf("%-4s %1d\n", P.taxon[i], keep[i]);
  //printf("\n");

  for (j = 0; j < ntre; j++)
    {
      PhylogenyAttachShuffle(P, S, attach);
      for (i = 0; i < P.nnodes; i++) keep[i] = 0;
      for (i = 0; i < S.srec[0]; i++) keep[P.t2n[attach[S.id[0][i]]]] = 1;

      Fy2newRec(Prune(P, keep));
    }

}

void SamplePrune(struct phylo P, sample S)
{
  int i, j;
  int *attach;
  int *keep;

  keep = ivector(0,P.nnodes-1);

  attach = ivector(0, S.ntaxa-1);
  AttachSampleToPhylo(S, P, attach);

  for (j = 0; j < S.nsamples; j++)
    {
      for (i = 0; i < P.nnodes; i++) keep[i] = 0;
      for (i = 0; i < S.srec[j]; i++) keep[P.t2n[attach[S.id[j][i]]]] = 1;

      Fy2newRec(Prune(P, keep));
    }

}

int CleanPhy(phylo P)
{
  // Modifies a phylo structure by relinking nodes `through' empty nodes
  // test with (((term1)empty1,a),(b,c,((term2)subempty2)empty2))root;

  // Also cuts the tail off.

  int i, at, tmp, mark;
  Droptail = 1;

  P = SetNodePointers(P);

  for (i = 0; i < P.nnodes; i++)
    {

      // if a one-daughter node, not yet excised:
      // not terminal && still with up && no rsister to ldown
      if ((P.ldown[i] != -99) && (P.up[i] != -99) && \
          (P.rsister[P.ldown[i]] == -99))
        {
          // printf("fixing %d\n", i);
          at = i;
    
          while (P.rsister[P.ldown[at]] == -99)
            {
              // printf("  doing %d\n", at);

              // fix branch lenghts
              P.bl[P.ldown[at]] += P.bl[at];
              // fix sisters
              // if left branch
              if (P.rsister[at] != -99) 
                {
                  P.rsister[P.ldown[at]] = P.rsister[at];
                }
              // if right branch 
              else if (P.rsister[P.ldown[P.up[at]]] != -99) 
                {
                  // find the dangling sister
                  mark = P.ldown[P.up[at]];
                  while (P.rsister[mark] != at)
                    {
                      mark = P.rsister[mark];
                    }
                  P.rsister[mark] = P.ldown[at];
                }
              // or if with another empty node above
              else if ( (P.rsister[at] == -99) && 
                        (P.rsister[P.ldown[P.up[at]]] == -99))
                {
                  // nothing needed
                }
              // fix daughter -> parent
              P.up[P.ldown[at]] = P.up[at];
              // fix parent -> daughter
              if (P.ldown[P.up[at]] == at)
                {
                  P.ldown[P.up[at]] = P.ldown[at];
                }

              tmp = at;
              // move point
              at = P.up[at];
              // erase old node
              P.up[tmp] = -99; 
            }
        }
    }      

  //for (i = 0; i < P.nnodes; i++)
  //  {
  //    printf("i %d, up %d, ldown %d, rsist %d, bl %f, name %s\n", 
  //         i, P.up[i], P.ldown[i], P.rsister[i], P.bl[i], P.taxon[i]);
  //  }

  if (FYOUT) FyOut(P) ;
  else Fy2newRec(P) ;

  return 1;
}
