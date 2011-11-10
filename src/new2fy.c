// new2fy.c Reads newick tree files

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

phylo New2fy(char filename[50]) {

  phylo Xtree;

  char *line; // array of characters from input line
  int lng = 0;
  char tmp[15];

  int nodei = -1;
  int leni = -1;
  char taxa[MAXTAXONLENGTH];
  int i, done, q, p;
  int count = 0; /* number of characters seen */
  int atn = 0;
  char iname[MAXTAXONLENGTH];
  int branchlengths = 0;

  /* character or EOF flag from input */
  int ch;
  int lbrack = 0;
  int rbrack = 0;
  int comma = 0;

  // Note - because the file parsed as a long string, not need
  // to add a line ending test.

  // PRE-READ
  Fn = fopen(filename, "r");
  if (Fn == NULL) {
    printf("Cannot open phylogeny file\n");
    exit(8);
  }

  while (1) {
    ch = fgetc(Fn);
    if (ch == 40)
      lbrack++;
    if (ch == 41)
      rbrack++;
    if (ch == 44)
      comma++;
    if (ch == EOF)
      break;
    ++count;
  }

  // printf("Number of characters in phylo.new is %d\n", count);
  // printf("Number of lbracks in phylo.new is %d\n", lbrack);
  // printf("Number of commas in phylo.new is %d\n", comma);
  // printf("Number of nodes/lines in phylo.fy should be %d\n", lbrack + comma + 1);

  fclose(Fn);

  if (lbrack != rbrack) {
    printf("Imbalanced parentheses in phylogeny file.\n");
    exit(EXIT_FAILURE);
  }

  // INITIALIZE
  Xtree.arenotes = 0;
  // Set the number of nodes in the tree
  Xtree.nnodes = lbrack + comma + 1;
  // Allocate size to structure:      
  Xtree.up = ivector(0, Xtree.nnodes-1);
  Xtree.down = imatrix(0, Xtree.nnodes-1, 0, comma); // set max poly = comma
  Xtree.noat = ivector(0, Xtree.nnodes-1);
  Xtree.depth = ivector(0, Xtree.nnodes-1);
  Xtree.bl = vector(0, Xtree.nnodes-1);
  Xtree.age = vector(0, Xtree.nnodes-1);
  Xtree.taxon = cmatrix(0, Xtree.nnodes-1, 0, MAXTAXONLENGTH);
  Xtree.ntaxa = 0;
  Xtree.termtaxa = 0;
  Xtree.taxalist = cmatrix(0, Xtree.nnodes-1, 0, MAXTAXONLENGTH);
  Xtree.t2n = ivector(0, Xtree.nnodes-1); // Maximum possible if all nodes are named

  for (p = 0; p < Xtree.nnodes; p++) {
    Xtree.noat[p] = 0;
    Xtree.depth[p] = 0;
    for (q = 0; q <= MAXTAXONLENGTH; q++)
      Xtree.taxon[p][q] = 0;
  }

  // Set size of line
  line = cvector(0, count); // number of boxes is count + 1 for final NULL


  // READ PROPERLY
  if ((Fn = fopen(filename, "r")) == NULL) {
    printf("Cannot open phylogeny file\n");
    exit(0);
  }

  fread(line, count, 1, Fn);
  fclose(Fn);
  line[count-1]=59;
  line[count]=0;

  // get length of newick - diff from count - there may be spaces or \n\n
  while (line[lng] != 59) {
    lng++;
  }

  // MOVE THROUGHT THE NEWICK FORMAT TREE CHARACTER BY CHARACTER
  i = 0;
  while (i < lng) {
    done = 0;

    // descend a branch, create lookup
    if (line[i] == 40) // "(" 
      {
        nodei++;

        if (leni == -1) {
          Xtree.depth[nodei] = 0;
        } else {
          Xtree.depth[nodei] = Xtree.depth[leni] + 1;
          Xtree.noat[leni]++;
          Xtree.down[leni][Xtree.noat[leni]-1] = nodei;
        }

        Xtree.up[nodei] = leni;
        strcat(Xtree.taxon[nodei], ".");
        leni = nodei;
        done = 1;
        i++;
      }

    // sibling taxa
    else if (strncmp(&line[i], ",", 1)==0) // ","
      {
        done = 1;
        i++;
      }

    // back up a node to len, keep track of locn with atn
    else if (strncmp(&line[i], ")", 1)==0) // ")"
      {
        leni = Xtree.up[leni];
        atn = Xtree.up[atn];
        done = 1;
        i++;
      }

    // Interior name  A-Za-z\-\_ and last one was )
    else if ( ( ( (line[i] >= 65) && (line[i] <= 90)  ) || \
                ( (line[i] >= 97) && (line[i] <= 122) ) || \
                (line[i] == 45)                         || \
                (line[i] == 95)                         ) && \
              ( strncmp(&line[i-1], ")", 1) == 0)          ) {
      // clear iname
      for (q = 0; q < MAXTAXONLENGTH; q++)
        iname[q] = 0;

      while ( (strncmp(&line[i], ":", 1) != 0) && \
              (strncmp(&line[i], ",", 1) != 0) && \
              (strncmp(&line[i], ")", 1) != 0) && \
              (strncmp(&line[i], "[", 1) != 0) && \
              (strncmp(&line[i], ";", 1) != 0) && \
              (strncmp(&line[i], "]", 1) != 0) ) 
        {
          strncat(iname, &line[i], 1);
          i++;
        }
      
      strcpy(Xtree.taxon[atn], iname);
      done = 1;
    }

    // NOTE - IGNORE IT
    else if (strncmp(&line[i], "[", 1) == 0) {
      while (strncmp(&line[i], "]", 1) !=0) {
        i++;
      }
      i ++; // to move over last one
      done = 1;
    }

    // branch length coming up
    else if (strncmp( &line[i], ":", 1) == 0) {

      char bl[30] = { 0 };

      i++;

      // 0-9 or .
      while ( ( (line[i] >= 48) && (line[i] <= 57)) || (line[i] == 46)) {
        strncat(bl, &line[i], 1);
        //printf("%s\n",bl); 
        i++;
      }

      if (bl[0] != 0)
        branchlengths = 1;
      //printf("%s\n",bl);  
      sscanf(bl, "%f", &Xtree.bl[atn]);
      //printf("%f is the bl of %d\n", Xtree.bl[atn], atn);
      done = 1;
    }

    // skip over whitespace before continuing to parse
    while ( (line[i] == ' ') || (line[i] == '\t') || (line[i] == '\n') || (line[i] == '\r')) {
      i++;
    }
                
    // default - it's a new taxon name
    // TODO - fix this. anytime you find a non ,:;()[] character, check to 
    // see if previous VALID character is , or ). If it is, or if character
    // is ', begin a new taxon name. If it's not, ignore. This should handle
    // whitespace, newlines and other garbage (but non-fatal) characters 
    // that currently crash the input.
    if (done == 0) {
      for (q = 0; q < 100; q++)
        taxa[q] = 0;
      // printf("i: %d   line[i]: %c   taxa: %s\n", i, line[i], taxa);
      strncat (taxa, &line[i], 1);
      // printf("i: %d   line[i]: %c   taxa: %s\n", i, line[i], taxa);
      i++;

      // KEEP ADDING MORE
      while ( (strncmp( &line[i], ",", 1) != 0) && \
              (strncmp( &line[i], ")", 1) != 0) && \
              (strncmp( &line[i], ":", 1) != 0) && \
              (strncmp( &line[i], "[", 1) != 0) ) 
        {
          strncat(taxa, &line[i], 1);
          i++;
        }

      // A NEW NAME MEANS A NEW NODE
      nodei++;
      atn = nodei;
      Xtree.depth[nodei] = Xtree.depth[leni] + 1;
      Xtree.up[nodei] = leni;
      strcpy(Xtree.taxon[nodei], taxa );
      Xtree.noat[leni]++;
      Xtree.down[leni][Xtree.noat[leni]-1] = nodei;
      // printf("%d ",nodei);

    }
  }

  // tidy up
  for (i = 0; i <= nodei; i++) {
    if (branchlengths == 0)
      Xtree.bl[i] = 1.0;
  }

  //    for (i = 0 ; i < Xtree.nnodes; i++) */
  //     { */
  //      printf("%d\t%d\t%d\t%d,%d,%d...\t%d\t%f\t%s\n",\ */
  //         i, Xtree.up[i], Xtree.noat[i], Xtree.down[i][0], Xtree.down[i][1],Xtree.down[i][2],Xtree.depth[i], Xtree.bl[i], Xtree.taxon[i]); */
  //     } */

  // build phylo metadata
  Xtree.ntaxa = 0;
  Xtree.termtaxa = 0;
  strcpy(Xtree.phyname, "phylo_from_file_");
  strcat(Xtree.phyname, filename);
  for (i = 0; i < Xtree.nnodes; i++) {
    if (strcmp(Xtree.taxon[i], ".") == 0) // now nothing is blank
      {
        if (MAKENODENAMES == 1) {
          // create default node names if so desired
          sprintf(tmp, "node%d", i);
          strcpy(Xtree.taxon[i], tmp);
        } else {
          strcpy(Xtree.taxon[i], "");
        }
      }

    //if node is named, increment count of node names
    if (strlen(Xtree.taxon[i]) > 0) {
      Xtree.ntaxa++;
    }
                
    //if this node is terminal, add taxon to taxalist
    if (Xtree.noat[i] == 0) {
      Xtree.termtaxa++;
      // create index of terminal taxa names .taxalist[terminal taxon #]
      strcpy(Xtree.taxalist[Xtree.termtaxa-1], Xtree.taxon[i]);
      // lookup to node from termtaxa#
      Xtree.t2n[Xtree.termtaxa-1] = i; // node # corresponding to term taxa #
    }
  }

  return Xtree;

}

