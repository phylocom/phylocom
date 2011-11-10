#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

void Comnode(phylo t1, phylo t2)
{
  int i, j, xnode, depth, common, matches1, matches2, p, q;
  int matchcode = 0;
  int t1md = 0; 
  int t2md = 0;
  int *active_n;
  int *matched1;
  int *matched2;
  int **comlist1; // contents: node number of a term taxon (use to lookup name)
  // indexed by node, and then an index 0...
  int *comlist1n; // the number of items in *comlist
  int **comlist2; // contents: node number of a term taxon (use to lookup name)
  // indexed by node, and then an index 0...
  int *comlist2n; // the number of items in *comlist
  char tmp[50];

  phylo Out[1];

  strcpy(t1.phyname, "Tree1");
  // t1.arenotes = 1;
  // t1.notes = cmatrix(0, t1.nnodes-1, 0, 49);
  strcpy(t2.phyname, "Tree2");
  // t2.arenotes = 1;
  // t2.notes = cmatrix(0, t2.nnodes-1, 0, 49);

  comlist1 = imatrix(0, t1.nnodes-1, 0, t1.termtaxa);
  comlist1n = ivector(0, t1.nnodes-1);
  comlist2 = imatrix(0, t2.nnodes-1, 0, t2.termtaxa);
  comlist2n = ivector(0, t2.nnodes-1);
  matched1 = ivector(0, t1.termtaxa-1);
  matched2 = ivector(0, t2.termtaxa-1);

  // printf("t1: nnodes: %d taxa1: %s\n", t1.nnodes, t1.taxalist[0]);
  // printf("t2: nnodes: %d taxa1: %s\n", t2.nnodes, t2.taxalist[0]);

  // make a full list of term taxa from each node:
  // tree1:
  for (i = 0; i < t1.nnodes; i++) 
    {
      if (t1.depth[i] > t1md) t1md = t1.depth[i];

      // need to clear the node names - they could cause confusion
      // when bladjing, beacause the default names are the same for both
      // trees
      if (t1.noat[i] > 0) strcpy(t1.taxon[i], ".");
      // these will get unmaned: see fy2new
    }
  active_n = ivector(0, t1.nnodes-1);

  for (i = 0; i < t1.nnodes; i++)
    {
      active_n[i] = 1;
      comlist1n[i] = 0; 

      for (depth = t1.depth[i]; depth <= t1md; depth++)
        {
          for (xnode = 0; xnode < t1.nnodes; xnode++)
            {
              if ( (t1.depth[xnode] == depth) && (active_n[xnode] == 1) )
                {
                  // first we check to see if we have reached a tip
                  if (t1.noat[xnode] == 0)
                    {
                      // test to see that it is a common spp to tree 2
                      common = 0;
                      for (j = 0; j < t2.termtaxa; j++)
                        {             
                          if (strcmp(t1.taxon[xnode], t2.taxalist[j]) == 0)\
                            common = 1;
                        }
                      if (common ==1)
                        { 
                          comlist1[i][comlist1n[i]] = xnode; 
                          comlist1n[i]++;
                          //printf("1 dependent to node %d is %s\n", 
                          //    i, t1.taxon[comlist1[i][comlist1n[i]-1]]);
                        }
                      active_n[xnode] = 0;
                    }
          
                  else
                    {
                      for (j = 0; j < t1.noat[xnode]; j++)
                        {
                          active_n[t1.down[xnode][j]] = 1;
                        }
                      active_n[xnode] = 0;
                    }
                }
            }
        }
    }

  // make a full list of term taxa from each node:
  // tree2:
  for (i = 0; i < t2.nnodes; i++) 
    {
      if (t2.depth[i] > t2md) t2md = t2.depth[i];
      // need to clear the node names - they could cause confusion
      // when bladjing, beacause the default names are the same for both
      // trees
      if (t2.noat[i] > 0) strcpy(t2.taxon[i], ".");
      // these will get unmaned: see fy2new
    }
  active_n = ivector(0, t2.nnodes-1);

  for (i = 0; i < t2.nnodes; i++)
    {
      active_n[i] = 1;
      comlist2n[i] = 0; 

      for (depth = t2.depth[i]; depth <= t2md; depth++)
        {
          for (xnode = 0; xnode < t2.nnodes; xnode++)
            {
              if ( (t2.depth[xnode] == depth) && (active_n[xnode] == 1) )
                {
                  // first we check to see if we have reached a tip
                  if (t2.noat[xnode] == 0)
                    {
                      // test to see that it is a common spp to tree 2
                      common = 0;
                      for (j = 0; j < t1.termtaxa; j++)
                        {             
                          if (strcmp(t2.taxon[xnode], t1.taxalist[j]) == 0)\
                            common = 1;
                        }
                      if (common ==1)
                        { 
                          comlist2[i][comlist2n[i]] = xnode; 
                          comlist2n[i]++;
                          //printf("2 dependent to node %d is %s\n", 
                          //      i, t2.taxon[comlist2[i][comlist2n[i]-1]]);
                        }
                      active_n[xnode] = 0;
                    }
          
                  else
                    {
                      for (j = 0; j < t2.noat[xnode]; j++)
                        {
                          active_n[t2.down[xnode][j]] = 1;
                        }
                      active_n[xnode] = 0;
                    }
                }
            }
        }
    }

  // Now, compare the nodes (once only for each pair)
  for (i = 0; i < t1.nnodes; i++)
    {
      for (j = 0; j < t2.nnodes; j++)
        {
          // initialize:
          for (p = 0; p < comlist1n[i]; p++) matched1[p] = 0;
          for (q = 0; q < comlist2n[j]; q++) matched2[q] = 0;

          // here's the crux - each one must have the same sublist to be valid
          for (p = 0; p < comlist1n[i]; p++)
            {
              for(q= 0; q < comlist2n[j]; q++)
                {
                  if (strcmp(t1.taxon[comlist1[i][p]], 
                             t2.taxon[comlist2[j][q]] ) ==0 )
                    {
                      matched1[p] = 1;
                      matched2[q] = 1;
                    }
                } 
            }
          // check matches
          matches1 = 0;
          matches2 = 0;
          for (p = 0; p < comlist1n[i]; p++)
            {
              if (matched1[p] == 1) matches1++;
            }
          for (q = 0; q < comlist2n[j]; q++)
            {
              if (matched2[q] == 1) matches2++;
            }
          if ( (matches1 == comlist1n[i]) && (matches2 == comlist2n[j]) && \
               (matches1 > 1) && (matches2 > 1))
            {
              sprintf(tmp, "match%d", matchcode++);
              //printf("tree1node%d matches with tree2node%d\n", i, j);

              // more correct, but need to just use nodenames for r8s:
              // strcat(t1.notes[i], tmp); strcat(t1.notes[i], " ");
              // strcat(t2.notes[j], tmp); strcat(t2.notes[j], " ");

              // for bladj (note the last one overrides the previous - this
              // will general be correct
              strcpy(t1.taxon[i], tmp);
              strcpy(t2.taxon[j], tmp);
            }
        }
    }

  // Unname the "." node names (a hack)
  for (i = 0; i < t1.nnodes; i++) 
    {
      if (!strcmp(t1.taxon[i], ".")) strcpy(t1.taxon[i], "");
    }
  for (i = 0; i < t2.nnodes; i++) 
    {
      if (!strcmp(t2.taxon[i], ".")) strcpy(t2.taxon[i], "");
    }
    

  printf("#NEXUS\n\nBEGIN TREES;\nTREE tree1 = ");
  Fy2newRec(t1);

  printf("TREE tree2 = ");
  Fy2newRec(t2);
  printf("END;\n");

  Out[0] = t1;
  //WriteNexus(Out, 1, InS, 0 , InC, 0 );

  Out[0] = t2;
  //WriteNexus(Out, 1, InS, 0 , InC, 0 );

}
