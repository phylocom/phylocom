#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "phylocom.h"
#include "nrutil.h"

void Bladj(phylo Intree)
{

  int i, j, q, z, l = 0;
  int *action;
  int *AgeFixed;
  char nameI[50];
  float ageI;
  int matched;
  char line[201];  // array of characters from input line
  int lineending;

  // Dimension things:
  action = ivector(0, Intree.nnodes -1);
  AgeFixed = ivector(0, Intree.nnodes -1);

  // Fix ages for nodes of terminal taxa in Phylogeny:
  for (i = 0; i < Intree.nnodes; i++)
    {
      // terminal nodes
      if (Intree.noat[i] == 0)
        {
          Intree.age[i] = 0.0;
          AgeFixed[i] = 1;
        }
      else
        {
          Intree.age[i] = 99999.9;
          AgeFixed[i] = 0;
        }
    }

  // pre-read
  lineending = whatnewline(INFILEA);
  
  //Fix node ages for nodes found in ages file
  Fa = fopen(INFILEA, "r");
  while (myfgets(line, 200, Fa, lineending) != NULL)
    {
      sscanf(line, "%s %f", nameI, &ageI); // string
      matched = 0;
      for (z = 0; z < Intree.nnodes; z++)
        {
          if (strcmp(Intree.taxon[z], nameI) == 0)
            {
              Intree.age[z] = ageI;
              AgeFixed[z] = 1;
              matched = 1;
            }
        }
    }
  fclose(Fa);

  //TODO will crash if no name/age for deepest node, need to add check for this

  // The algorithm:
  // 1. create network of fixed age nodes between the root and the other
  //    fixed age nodes, choosing the order of nodes to operate on first
  //    using age, then number of intervening nodes.

  for (i = 0; i < Intree.nnodes; i++)
    {
      q = 0;
      for (j = i+1; j < Intree.nnodes; j++)
        {
          // find all the line-of-site ages
          action[q] = 0;
      
          // correct for errors in ages
          if ((LineOfSight(Intree, AgeFixed, i, j) == 1) && \
              (Intree.age[j] >= Intree.age[i])) AgeFixed[j] = 0;
      
          if ((AgeFixed[j] == 1) && (LineOfSight(Intree, AgeFixed, i, j) == 1))
            {
              // printf("%d+%d  ",i, j);
              action[q] = j;
              q++;
            }
        }
      
      // Now sort the action  
      SortAction(Intree, action, q, i);
      
      // Adjust lengths
      for (l = 0; l < q; l++)
        {
          Adjust(Intree, AgeFixed, i, action[l]);
          //printf("i%d action%d\n",i, action[l]);
        }
    }
  Fy2newRec(Intree);
}

//void Polytom(phylo Intree)
//{
//  // need to: work out which nodes (n) are polytomies
//  // workout what the % shortening of dependent nodes should be
//  // find the shortest dependent segment and shorten by amount x
//  // shorten all others by absolute amount x
//  // extend the bl(n)
//
//  int i;
//
//  for (i = 0; i <= MaxNode; i++)
//  {
//    if (Intree.noat[Intree.up[i]] > 2)
//      {
//        // Adjust the BL to reflect the mean
//      }
//  }
//}

void Adjust(phylo Intree, int AgeFixed[], int inner, int outer)
{
  // minor untidiness: the BLs are reassigned (to same value) numerous times
  int z, y, x, segs;
  float segBL;

  // find the nearest proximal fixed age
  z = Intree.up[outer];
  while(AgeFixed[z] != 1)
    {
      //TODO this could cause error - need to make sure deepest node has name/age
      z = Intree.up[z];
    }
  //printf("action=%d inner=%d\n", outer, z);

  // Hwo many segments between outer and z?
  segs = Intree.depth[outer] - Intree.depth[z];
  segBL = (Intree.age[z] - Intree.age[outer]) / (float) segs;

  // printf("   segs=%d segBL=%f\n", segs, segBL);

  // Apply results
  y = outer;
  Intree.bl[y] = segBL;
  y = Intree.up[y]; x = 1; 
  while(y != z)
    {
      Intree.bl[y] = segBL;
      Intree.age[y] = Intree.age[outer] + (segBL * (float) x);
      AgeFixed[y] = 1;
      y = Intree.up[y]; x++;
    }
}

void SortAction(phylo Intree, int Array[], int Nia, int start)
{
  int i;
  float swap;
  int aSwap = 1;

  while (aSwap == 1)
    {
      aSwap = 0;

      for(i = 0; i < Nia -1; i++) /* note up to one less than final Total... */
        {
          if (Intree.age[Array[i]] < Intree.age[Array[i+1]])
            {
              swap = Array[i+1];
              Array[i+1] = Array[i];
              Array[i] = swap;
              aSwap = 1;
            }
          else if ((Intree.age[Array[i]] == Intree.age[Array[i+1]]) && \
                   ((Intree.depth[Array[i]] - Intree.depth[start]) < (Intree.depth[Array[i+1]] - Intree.depth[start])))
            {
              swap = Array[i+1];
              Array[i+1] = Array[i];
              Array[i] = swap;
              aSwap = 1;
            }
        }
    }
}

int LineOfSight(phylo Intree, int AgeFixed[], int a, int b)
{
  int inner, outer;
  int z = 1;

  if (Intree.depth[a] > Intree.depth[b]) 
    {
      inner = b; outer = a;
    }
  else
    {
      inner = a; outer = b;
    }

  while ((inner != outer) && (outer != -1))  
    {
      outer = Intree.up[outer]; 
      if ((AgeFixed[outer] == 1) && (outer != inner)) z = 0;
    }

  if (outer == -1) z = 0;
  return z;
}

