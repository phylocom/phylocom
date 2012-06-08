/*  phylomatic.c  */

// This program is free software; you can redistribute it and/or
// modify it under the terms of the BSD 2-Clause License
//
// http://www.opensource.org/licenses/bsd-license.php
//
// Copyright (c) 2005-2011, Campbell Webb <cwebb@oeb.harvard.edu>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.


// TODO: create a phylogenetic hierarchy from the nested clade names in TAXA
//       even if it does not exist in PHYLO - need to loop the search through
//       the new nodes.  Also need a switch to stop (singletaxon)name for
//       apps that cant take it - ahh - no, the pruning alg fixes this

/* INCLUDE HEADERS ------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "phylocom.h"
#include "nrutil.h"

int LOWERCASETAXA = 0;
int Droptail = 1;  // no root tails from phylomatic
int FYOUT = 0;

struct taxa ReadTaxa();

typedef struct taxa {
  char ***str; // nested taxa strings
  int ntx; // number of taxa
  int maxnest; // max number of nestedness
} taxa;


int main(int argc, char *argv[])
{
  struct phylo M;
  struct phylo P;
  struct taxa T;
  int i, j, k, p, next, last, nexttl;
  int nmatched = 0;
  int *matched;
  int *keep;
  int argx;
  char PHYLOFILE[100];
  char TAXAFILE[100];
  float totbl;

  MAKENODENAMES = 0;
  strcpy(TAXAFILE, "taxa");
  strcpy(PHYLOFILE, "phylo");
  Debug = 0;

  // set the buffer to zero, so output is written directly to file:
  setbuf(stdout, NULL);  

  // pick up the global switches:
  for (argx = 0; argx <  argc; argx++)
    {
      if (strcmp(argv[argx], "-f") == 0)
        {
          sscanf(argv[argx+1], "%s", PHYLOFILE);
        }
      if (strcmp(argv[argx],"-t") == 0)
        {
          sscanf(argv[argx+1], "%s", TAXAFILE);
        }
      if (strcmp(argv[argx],"-n") == 0)
        {
          MAKENODENAMES = 1;
        }
      if (strcmp(argv[argx],"-l") == 0)
        {
          LOWERCASETAXA = 1;
        }
      if (strcmp(argv[argx],"-d") == 0)
        {
          Debug = 1;
        }
      if (strcmp(argv[argx],"-y") == 0)
        {
          FYOUT = 1;
        }
      if ((strcmp(argv[argx], "-h") == 0) || (strcmp(argv[argx], "--help") == 0))
        {
          printf("\n    ======================================================\n\n");
          printf("                     P H Y L O M A T I C\n");
          printf("                Cam Webb <cwebb@oeb.harvard.edu>\n");
          printf("\n    ======================================================\n\n");
          printf("Usage: phylomatic [-t -f -n -h -l]\n\n");
          printf("Switches:\n");
          printf("  -t FILENAME use this as the taxa file\n");
          printf("  -f FILENAME use this as the phylo file\n");
          printf("  -n    label all nodes with default names\n");
          printf("  -h    Help information\n");
          printf("  -l    Convert all chars in taxa file to lowercase\n\n");
          return 1;
        }
    }


  M = ReadPhylogeny(PHYLOFILE);
  T = ReadTaxa(TAXAFILE);

  M.age = vector(0, M.nnodes-1);
  AgeNodes(M);


  // create the new phylo, P:
  P.nnodes = M.nnodes + (T.ntx * T.maxnest); // to allow for lots of new hiera
  P.ntaxa = P.nnodes;
  P.up = ivector(0, P.nnodes-1);
  P.depth = ivector(0, P.nnodes-1);
  P.bl = vector(0, P.nnodes-1);
  P.taxon = cmatrix(0, P.nnodes-1, 0, MAXTAXONLENGTH);
  P.age = vector(0, P.nnodes-1);
  P.noat = ivector(0, P.nnodes-1);
  P.taxalist = cmatrix(0, P.ntaxa-1, 0, MAXTAXONLENGTH);
  P.t2n = ivector(0, P.ntaxa-1);

  for (i = 0; i < M.nnodes; i++)
    {
      P.up[i] = M.up[i];
      P.noat[i] = M.noat[i];
      P.depth[i] = M.depth[i];
      P.bl[i] = M.bl[i];
      P.age[i] = M.age[i];
      strcpy(P.taxon[i], M.taxon[i]);
    }
  for (i = 0; i < M.ntaxa; i++)
    {
      strcpy(P.taxalist[i], M.taxalist[i]);
      P.t2n[i] = M.t2n[i];
    }
  
  matched = ivector(0, T.ntx-1); for (i=0;i<T.ntx;i++) matched[i]=0;
  keep    = ivector(0, P.nnodes -1); 
  for (i=0;i<P.nnodes;i++) keep[i]=0;
  next = M.nnodes;
  nexttl = M.ntaxa;

  // - Main phylomatic algorithm
  // -- for each taxon in the `taxa` file
  for (i = 0; i < T.ntx; i++)
    {
      // --- climb through the levels of nesting from terminal to basal
      for (j = 0; j < T.maxnest; j++)
        {
          if (matched[i] != 1)
            {
              for (k = 0; k < next; k++) // looping through all nodes
                {
                  // ---- compare with the taxon name at every node in the 
                  // ---- megatree
                  if (strcmp(T.str[i][j], P.taxon[k]) == 0)
                    {
                      // printf("MATCH t%d+n%d:p%d %s\n", i, j, k, M.taxon[k]);
                      // ----- if matches, add to megatree:
                      // ------ if terminal matches terminal, join new taxon
                      // ------ to parent node of matching taxon in megatree
                      // terminal to terminal - perfect
                      if ((j == 0) && (P.noat[k] == 0))
                        {
                          P.bl[next] = P.bl[k];
                          P.up[next] = P.up[k];
                          P.depth[next] = P.depth[k];
                          P.noat[P.up[k]]++;
                          strcpy(P.taxon[next], T.str[i][0]);
                          strcpy(P.taxalist[nexttl], T.str[i][0]);
                          keep[next] = 1;
                          P.t2n[nexttl] = next;
                          P.age[next] = 0.0;
                          next++; nexttl++;
                          matched[i] = 1; nmatched++;
                        }
                      // interior to terminal
                      // ------ if interior matches interior, join new node path
                      // ------ to parent node of matching taxon in megatree

                      else if ((j > 0) && (P.noat[k] == 0))
                        {
                          last = k;
                          totbl = P.bl[last];
                          // indent terminals
                          P.bl[k] = totbl / (float) (j+1);
                          P.age[k] = P.age[P.up[k]] - P.bl[k];
                          for (p = j-1; p >= 0; p--)
                            {
                              P.bl[next] = totbl / (float)(j+1);
                              P.age[next] = P.age[k] - \
                                (totbl / (float) (j+1));
                              P.up[next] = last;
                              P.depth[next] = P.depth[last]+1;
                              P.noat[last]++;
                              strcpy(P.taxon[next], T.str[i][p]);
                              strcpy(P.taxalist[nexttl], T.str[i][p]);
                              if (p == 0) keep[next] = 1;
                              P.t2n[nexttl] = next;
                              last = next;
                              next++; nexttl++;
                            }
                          matched[i] = 1; nmatched++;
                        }
                      // interior to multi interior
                      else
                        {
                          last = k;
                          totbl = P.age[last];
                          // n added nodes = j
                          for (p = j-1; p >= 0; p--)
                            {
                              // make the interior first
                              P.bl[next] = totbl / (float) j;
                              P.age[next] = P.age[last] - \
                                (totbl / (float) j);
                              P.up[next] = last;
                              P.depth[next] = P.depth[last]+1;
                              P.noat[last]++;
                              strcpy(P.taxon[next], T.str[i][p]);
                              strcpy(P.taxalist[nexttl], T.str[i][p]);
                              if (p == 0) keep[next] = 1;
                              P.t2n[nexttl] = next;
                              last = next;
                              next++; nexttl++;

                            }
                          matched[i] = 1; nmatched++;
                        }
                      break;
                    }
                }
            } 
        }
    }
  
  if (Debug)
    {
      for (i = 0; i < P.nnodes; i++)
        {
          printf("%4d  %4d  %4d  %d  %5.1f  %5.1f  %-10s  %2d\n", i, P.up[i], P.noat[i], P.depth[i], P.age[i], P.bl[i],P.taxon[i], keep[i]);
        }
    }

  // this is it:

  if (FYOUT)
    {
      for (i = 0 ; i < P.nnodes; i++)
        {
          printf("%d\t%d\t%d\t\t%d\t%f\t%s\n",  \
                 i, P.up[i], P.noat[i],         \
                 P.depth[i],                    \
                 P.bl[i], P.taxon[i]);
        }
    }
  else
    {
      Fy2newRec(Prune(P, keep));
    }

  
  // Missing taxa
  if (nmatched < T.ntx)
    {
	  // fprintf(stderr, "\f");
      // fprintf(stderr, "comment=\"NOTE: %d taxa not matched\"\n", (T.ntx-nmatched));
	  fprintf(stderr, "NOTE: %d taxa not matched:\n", (T.ntx-nmatched));
    }
  P.nnodes = next;
  P.ntaxa = nexttl;

  // List non matches
  for (i = 0; i < T.ntx; i++)
    {
      if (!matched[i])
        {
          // fprintf(stderr, "extra_taxon=%s\n", T.str[i][0]);
		  fprintf(stderr, "%s\n", T.str[i][0]);
        }
    }

  free_c3d(T.str, 0, T.ntx-1, 0, T.maxnest-1, 0, MAXTAXONLENGTH);
  free_cmatrix(P.taxon, 0, P.nnodes-1, 0, MAXTAXONLENGTH);
  free_cmatrix(P.taxalist, 0, P.nnodes-1, 0, MAXTAXONLENGTH);

  return 1;
}

struct taxa ReadTaxa(char filename[100])
{
  FILE *taxafile;
  char line[1000];
  int slash, i, j, k;
  int maxslash  = 0;
  int linesread = 0;
  struct taxa T;
  int lineending;

  if ( (taxafile = fopen(filename, "r")) == NULL)
    {
      printf("Cannot find file: %s\n", filename);
      exit(EXIT_FAILURE);
    }

  lineending = whatnewline(filename);

  // preread to define slash number
  while (myfgets(line, 1000, taxafile, lineending) != NULL)
    {
      slash = 0;
      i = 0;
      while (line[i] != 0)
        {
          if (line[i] == 47) slash++;
          i++;
        }
      //printf("len %d: %s\n", (int) strlen(line), line);
      if (slash > maxslash) maxslash = slash;
      linesread++;
    }
  fclose(taxafile);
  //printf("slashes:%d lines:%d\n", maxslash, linesread);

  T.ntx = linesread;
  T.maxnest = maxslash+1;
  T.str = c3d(0, T.ntx-1, 0, T.maxnest-1, 0, MAXTAXONLENGTH);
  
  // initialize - important!
  for (i = 0; i < T.ntx; i++) {
    for (j =  0; j < T.maxnest; j++) {
      strcpy(T.str[i][j], ""); } }
 
  // now read:
  linesread = 0;
  taxafile = fopen(filename, "r");
  while (myfgets(line, 1000, taxafile, lineending) != NULL)
    // while (fgets(line, 1000, taxafile) != NULL)
    {
      i = (int) strlen(line) - 2; // not the newline

      // Convert to lowercase
      if (LOWERCASETAXA) 
        {
          for (k = 0; k < strlen(line); k++) line[k] = tolower(line[k]);
        }
      // printf("%s", line);

      j = i;
      k = 0;
      while (i >= -1) // newline
        {
          if ((line[i] == 47) || (i == -1) ) 
            {
              strncat(T.str[linesread][k], &line[i+1], j-i);
              i--;
              j = i;
              k++;
            }

          //if (line[i] == 47) slash++;
          i--;
        }
      linesread++;
    }
  fclose(taxafile);

  //  for (i = 0; i < T.ntx; i++) {
  //  for (j = 0; j < T.maxnest; j++) {
  //    printf("%d %d %s\n", i, j, T.str[i][j]); } }

  return T;

}


