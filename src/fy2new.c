// fy2new.c - Outputs a newick file with BLs and internal names.

//TODO add switch to allow internal names to be excluded (e.g. for CACTUS)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

void Fy2new(struct phylo O)
{
  int i, d, j, qq;
  char blpr[30]; // for bl printing
  //char **ClusterAt;
  // problem with huge matrices - this is a temp workaound
  //TODO fix this?
  char *ClusterAt[100000]; // big file!
  int allocated[100000];
  char **InteriorName;
  int HighestNodeNo;
  int MaxDepth = 0;
  char *tmp;

  InteriorName = cmatrix(0, O.nnodes-1, 0, MAXTAXONLENGTH);
  tmp = cvector(0, (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + MAXNOTELENGTH + 3)));

  for (i = 0; i < O.nnodes; i++)
    {
      allocated[i] = 0;

      // if (strcmp(O.notes[i], "") != 0) printf("%s\n", O.notes[i]);
      if ((strcmp(O.taxon[i], ".") != 0) && (O.noat[i] == 0))
        {   
          ClusterAt[i] = cvector(0, (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + MAXNOTELENGTH + 3))); 
          allocated[i] = 1; 
          // printf("%d allocated\n", i);
      
          strcpy(ClusterAt[i], O.taxon[i]);
          strcat(ClusterAt[i], ":" );
          sprintf(blpr, "%f", O.bl[i]);
          strcat(ClusterAt[i], blpr);
        }
      //      else
      //    {
      //      strcpy(ClusterAt[i], "" );
      //    }

      if ((strcmp(O.taxon[i], ".") != 0) && (O.noat[i] > 0))
        {
          strcpy(InteriorName[i], O.taxon[i]);
        }

      if (MaxDepth < O.depth[i]) MaxDepth = O.depth[i];
    
    }

  HighestNodeNo = O.nnodes - 1;

  //  for (i = 0; i < 8; i ++)
  //  {
  //    printf("%d\t%s\n", i, ClusterAt[i]);
  //  }

  // printf("-------\n");

  for (d = MaxDepth; d >= 0; d--)
    {
      for (i = 0; i<= HighestNodeNo; i++)
        {
          if (O.depth[i] == d)
            {
              // allocate space for clusterat
              if (allocated[i] == 0)
                {
                  ClusterAt[i] = cvector(0, (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + MAXNOTELENGTH + 3))); 
                  allocated[i] = 1; 
                  // printf("%d allocated (up=%d)\n", i, O.up[i]);
                  strcpy(ClusterAt[i], "");
                }
              if ((allocated[O.up[i]] == 0) && (O.up[i] != -1))
                {
                  ClusterAt[O.up[i]] = cvector(0, (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + MAXNOTELENGTH + 3))); 
                  allocated[O.up[i]] = 1; 
                  // printf("%d allocated from %d (switch = %d)\n", O.up[i], i,  allocated[O.up[i]]);
                  strcpy(ClusterAt[O.up[i]], "");
                }
          
              for (j = 0; j <= HighestNodeNo; j++)
                {
                  if (O.depth[j] == d)
                    {
                      // printf("%d\t%d\t%d\t%s\n", d, i, j, ClusterAt[j]);

                      // if there are two full terminal nodes with nothing
                      // above
                      if ((i != j) && \
                          (strcmp(ClusterAt[O.up[i]], "") == 0) && \
                          (strcmp(ClusterAt[i] , "") != 0) && \
                          (strcmp(ClusterAt[j], "" ) != 0) && \
                          (O.up[i] == O.up[j]))
                        {
                          strcpy(ClusterAt[O.up[i]], "(");
                          strcat(ClusterAt[O.up[i]], ClusterAt[i]);
                          strcat(ClusterAt[O.up[i]], ",");
                          strcat(ClusterAt[O.up[i]], ClusterAt[j]);
                          strcat(ClusterAt[O.up[i]], ")");
                          strcpy(ClusterAt[i] ,"");
                          strcpy(ClusterAt[j] ,"");

                          // new bit
                          if (strcmp(InteriorName[O.up[i]], "") != 0)
                            {
                              strcat(ClusterAt[O.up[i]], InteriorName[O.up[i]]);
                            }

                          strcpy(ClusterAt[O.up[i]], ClusterAt[O.up[i]]);
                          strcat(ClusterAt[O.up[i]], ":");
                          sprintf(blpr, "%f",  O.bl[O.up[i]] );
                          strcat(ClusterAt[O.up[i]], blpr); 
                          if ((O.arenotes > 0) && \
                              (strcmp(O.notes[O.up[i]], "") != 0))
                            {
                              strcat(ClusterAt[O.up[i]], "[%note = 'string:");
                              strcat(ClusterAt[O.up[i]], O.notes[O.up[i]]);
                              strcat(ClusterAt[O.up[i]], "']");
                              // if (strcmp(O.notes[i], "") != 0) printf("%s\n%s\n", O.notes[i], ClusterAt[O.up[i]]);
                            }
                        }

                      // if there is now a cluster above - adding to a poly
                      if ((i != j) && \
                          (strcmp(ClusterAt[O.up[i]], "") != 0) && \
                          (strcmp(ClusterAt[j], "" ) != 0) && \
                          (O.up[i] == O.up[j]))
                        {
                          // find the last )
                          for (qq = strlen(ClusterAt[O.up[i]]); qq > 0; qq--)
                            {
                              if (strncmp(&ClusterAt[O.up[i]][qq], ")", 1) \
                                  == 0) break;
                            }
                          // printf("%d\t%d\t%d\t%s\n", i, j, qq,  ClusterAt[O.up[i]]);
                          strcpy(tmp, ClusterAt[O.up[i]]);
                          // printf("> > > %s\n", tmp);
                          strcpy(ClusterAt[O.up[i]], "");
                          strncat(ClusterAt[O.up[i]], tmp, qq); // OW! changed from -1 to -0
                          strcat(ClusterAt[O.up[i]], ",");
                          strcat(ClusterAt[O.up[i]], ClusterAt[j]);
                          strcat(ClusterAt[O.up[i]], ")");
                          strcat(ClusterAt[O.up[i]], &tmp[qq+1]);
          
                          strcpy(ClusterAt[j], "");
 
                        }
                    }
                }
            }
        }

      //leftovers:
      for (i = 0; i<= HighestNodeNo; i++)
        {
          if ((O.depth[i] == d) && \
              (strcmp(ClusterAt[i], "") != 0) && (d != 0))
            {
              strcpy(ClusterAt[O.up[i]], ClusterAt[i]);
              strcpy(ClusterAt[i], "");
            }
          if ((O.depth[i] == d) && (i != 0))
            {
              // free memory
              if (allocated[i] == 1)
                {
                  free_cvector(ClusterAt[i], 0, (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + MAXNOTELENGTH + 3)));
                  // printf("               %d deallocated\n", i);
                }
          

            }
        }
    }

  //remove the root BL
  //for (qq = length(ClusterAt[0]); qq > 0; qq--)
  //  {
  //    if (substr(ClusterAt[0], qq, 1) == ")")
  //    {break};
  //  }

  printf("%s;\n", ClusterAt[0]);

  //  free_cmatrix(ClusterAt, 0, O.nnodes-1, 0, 
  //          (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + MAXNOTELENGTH + 3)));


  free_cmatrix(InteriorName, 0, O.nnodes-1, 0, MAXTAXONLENGTH);
  free_cvector(tmp, 0, (O.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + \
                                    MAXNOTELENGTH + 3)));

}

// NOW! The recursive way ;-)


struct phylo SetNodePointers(struct phylo A)
{
  int mark, x, y;
  int *first;

  // Dimension:

  A.ldown = ivector(0, A.nnodes-1);
  A.rsister = ivector(0, A.nnodes-1);
  first = ivector(0, A.nnodes-1);

  for (x = 0; x < A.nnodes; x++)
    {
      A.ldown[x] = -99;
      A.rsister[x] = -99;
      first[x] = 1;
    }

  for (x = 0; x < A.nnodes; x++)
    {
      // starting at terms
      if (A.noat[x] == 0)
        {
          y = x;
          // while not yet at the root
          while (y != 0)
            {
              // is this the first daughter in new structure?
              if (A.ldown[A.up[y]] == -99) 
                {
                  A.ldown[A.up[y]] = y;
                }
              // if not, find the dangling sister
              else
                {
                  // start at ldown
                  mark = A.ldown[A.up[y]];
                  // move to node with an empty rsister
                  while (A.rsister[mark] != -99)
                    {
                      mark = A.rsister[mark];
                    }
                  A.rsister[mark] = y;
                }

              // test for refollowing old routes
              if (first[A.up[y]] == 1) {y = A.up[y]; first[y] = 0;}
              else break;
              
            }
        }
    }

  //for (x = 0; x < A.nnodes; x++)
  //  {
  //    printf("%d\t%d\t%d\t%d\t%s\n", x, A.up[x], A.ldown[x], A.rsister[x],A.taxon[x]);
  //  }

  free_ivector(first, 0, A.nnodes-1);

  return A;
}


void Fy2newRec(struct phylo A)
{
  char *tmp;
  // printf(">>> %d\n\n", A.nnodes);
  tmp = cvector(0, (A.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + \
                                MAXNOTELENGTH + 3)));

  // add ldown, rsister structures
  A = SetNodePointers(A);

  printf("%s;\n", downPar(A, 0, tmp));
  free_cvector(tmp, 0, (A.nnodes * (MAXBLLENGTH + MAXTAXONLENGTH + \
                                    MAXNOTELENGTH + 3)));

}

char *downPar(struct phylo A, int atn, char *tmp)
{
  int x;
  int which = 0;
  char blstring[20];
  char **tmpnext;

  tmpnext = cmatrix(0, A.noat[atn]-1, 0, (A.nnodes * (MAXBLLENGTH + \
                                                      MAXTAXONLENGTH + MAXNOTELENGTH + 3)));

  if (A.noat[atn] == 0)
    {
      strcpy(tmp, A.taxon[atn]);
      if (!NoBL) 
        {
          strcat(tmp, ":");
          sprintf(blstring, "%f", A.bl[atn]);
          strcat(tmp, blstring);
        }
    }
  else
    {
      x = A.ldown[atn];
      strcpy(tmp, "(");
      strcat(tmp, downPar(A,x, tmpnext[which]));
    
      x = A.rsister[x]; which++;

      while (x != -99)
        {
          strcat(tmp, ",");
          strcat(tmp, downPar(A, x, tmpnext[which]));
          x = A.rsister[x]; which++;
        }
      strcat(tmp, ")");
      strcat(tmp, A.taxon[atn]);
      if (!NoBL) 
        {
          // root tail causes probs in some other apps 
          if ((atn != 0) || (Droptail == 0))
            {
              strcat(tmp, ":");
              sprintf(blstring, "%f", A.bl[atn]);
              strcat(tmp, blstring);
            }
        }
      // notes:
      if ((A.arenotes > 0) && (strcmp(A.notes[atn], "") != 0))
        {
          strcat(tmp, "[%note = 'string:");
          strcat(tmp, A.notes[atn]);
          strcat(tmp, "']");
        }
    }
  
  free_cmatrix(tmpnext, 0, A.noat[atn]-1, 0, (A.nnodes * (MAXBLLENGTH   \
                                                          + MAXTAXONLENGTH + MAXNOTELENGTH + 3)));

  return tmp;

}
