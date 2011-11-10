/*  ecomain.c  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "ecovolve.h"
#include "nrutil.h"

/* MAIN PROGRAM ---------------------------------------------------------- */

int main(int argc, char *argv[])
{
  int l, t, ll, argx,   tout, aliveout;
  int lroute;
  int alive = 0;
  float p_SPECIATE = P_SPECIATE;
  float p_EXTINCT  = P_EXTINCT;
  float newp_EXTINCT;
  //char CHANGE[11] = {'3', '2', '1', '1', '0', '0', '0', '0', '0', '0', '\0'};
  char CHANGE[11] = "3211000000";
  int stop = 0;
  char STOPPER[11] = "time";
  int LINEAGES = 32;
  int BAILT = 500;
  int chi;

  PRUNED = 0;
  Droptail = 0;
 
  setbuf(stdout, NULL);  
  // srandom(time(NULL));
  //srandom(2);
  srandom(getpid()); // need this because in a shell-script the program 
  // may run several times

  MAX_TIME = MAXTIME;
  OUT_MODE = OUTMODE;
  NCHAR = 5;

  // read the arguments
  for (argx = 0; argx <  argc; argx++)
    {
      // runs
      if (strcmp(argv[argx], "-s") == 0)
        {
          sscanf(argv[argx+1], "%f", &p_SPECIATE);
        }
      if (strcmp(argv[argx], "-e") == 0)
        {
          sscanf(argv[argx+1], "%f", &p_EXTINCT);
        }      
      if (strcmp(argv[argx], "-t") == 0)
        {
          sscanf(argv[argx+1], "%d", &MAX_TIME);
        } 
      if (strcmp(argv[argx], "-d") == 0)
        {
          TAPER = 1;
          sscanf(argv[argx+1], "%d", &TAPERFACT);
        }
      if (strcmp(argv[argx], "-m") == 0)
        {
          sscanf(argv[argx+1], "%d", &OUT_MODE);
        }
      if (strcmp(argv[argx], "-c") == 0)
        {
          sscanf(argv[argx+1], "%s", CHANGE);
        }
      if (strcmp(argv[argx], "-p") == 0)
        {
          PRUNED = 1;
        }
      if (strcmp(argv[argx], "-l") == 0)
        {
          strcpy(STOPPER, "lineages");
          sscanf(argv[argx+1], "%d", &LINEAGES);
        }

      if (strcmp(argv[argx], "-x") == 0)
        {
          COMPETE = 1;
        }

      if ((strcmp(argv[argx], "-h") == 0) || (strcmp(argv[argx], "--help") == 0))
        {
          printf("\n    ======================================================\n\n");
          printf("                       E C O V O L V E\n");
          printf("    PHYLOGENY AND TRAIT SIMULATOR WITH ECOLOGICAL PROCESSES\n");
          printf("                Cam Webb <campbell.webb@yale.edu>\n");
          printf("\n    ======================================================\n\n");
          printf("Usage: ecovolve [-s -e -t -m -c -p -l -d -h]\n\n");
          printf("Switches:\n");
          printf("  -s    prob of speciation per unit time[0.05]\n");
          printf("  -e    prob of extinction per unit time[0.01]\n");
          printf("  -t    time units to simulate over [100]\n");
          printf("  -m    output mode (2 = LTT; 3 = newick) [3]\n");
          printf("  -c    prob envelope for char change, string of 10 ints [3211000000]\n");
          printf("  -l    stop simulation after this number of extant lineages [NO]\n");
          printf("  -p    output phylogeny pruned only for extant taxa [NO]\n");
          printf("  -d    taper character change by exp(-time / INT )\n");
          printf("  -x    simulate competition, with trait proximity increasing extinction\n");
          printf("  -h    help\n\n");
          printf("View it!:\n");
          printf("  $ ecovolve > ecovolve.phylo\n");
          printf("  $ phylocom makenex -f ecovolve.phylo -t ecovolve.traits \\\n");
          printf("             -s ecovolve.sample > ecovolve.nex\n");
          printf("  (open and trace character in Mesquite with BLs proportional)\n\n");
          return 1;
        }
    }

  // initialize dynamic arrays
  CharMin = vector(0, NCHAR-1);
  CharMax = vector(0, NCHAR-1);
  CharMean = vector(0, NCHAR-1);
  CharSumXSqr = vector(0, NCHAR-1);
  CharCMin = vector(0, NCHAR-1);
  CharCMax = vector(0, NCHAR-1);
  CharCMean = vector(0, NCHAR-1);
  CharCSumXSqr = vector(0, NCHAR-1);
  Char_ltc = f3tensor(0, MAXNODES-1, 0, MAXTIME-1, 0, NCHAR-1);
  Char_nc = matrix(0, MAXNODES-1, 0, NCHAR-1);

  // initalize lineage/time space:
  for (l = 0; l < MAXNODES; l++)
    {
      for (t = 0; t < MAX_TIME+1; t ++)
        {
          Living_lt[l][t] = 0;
          Censor_lt[l][t] = 0;
        }
    }
    
  // initialize starting lineage
  Node = 0; // AHA A problem here!!!! switch to 0 must change later
  Lineage = 0;
  Active_l[Lineage] = 1;
  Living_lt[Lineage][0] = 1;
  NodeUp_l[Lineage] = Node;
  for (chi = 0; chi < NCHAR; chi++)
    {
      Char_ltc[Lineage][0][chi] = 0.0;
      Char_nc[Node][chi] = 0.0;
    }
  Time_n[Node] = 0;
  Name_n[Node] = 0;
  EcoDist_l[Lineage] = 10000;
  Name = 0;

  // Prepare change vector
  charch = MakeChange(CHANGE);

  // for each time slice...
  t = 0;
  while (stop == 0)
    {
      
      // Determine competitive effect on extinction (before new evolution
      // and speciation)  simplest model - distance to nearest phenotype
      Compete(t);
                        
      // store the last lineage number, so that it doesnt get confused by new species
      ll = Lineage;

      alive = 0;
      // for each active lineage
      for (l = 0; l <= ll; l++)
        {
          if (Active_l[l] == 1)
            {
              alive++;

              // evolve Character - NB gradual on lineage
              CharChange(l, t);

              // determine competetive environments
              newp_EXTINCT = p_EXTINCT;
              // vvv EDIT HERE FOR COMPETITION ALG vvv
              if (COMPETE)
                {
                  // the modified extinction probability is 
                  newp_EXTINCT = p_EXTINCT * (5.0 / ((float) EcoDist_l[l] \
                                                     + 1.0));
                  // printf("%f\n",newp_EXTINCT);
                }
              // ^^^                               ^^^

              if (OUT_MODE == 0)
                {
                  printf("      t %d\tl %d\t\t\t\t\tc %f\te %d\n",  \
                         t, l, Char_ltc[l][t][0], EcoDist_l[l]);
                }

              // check for speciation
              if ( ((float) random() / (float) (RAND_MAX+1.0)) < p_SPECIATE)
                {
                  Speciate(l, t);
                  alive++;
          
                }

              // check for extinction, if the lineage has not speciated,
              // and if the lineage is not the very first
              else if (( ((float) random() / (float) (RAND_MAX+1.0))      \
                         < newp_EXTINCT)                  \
                       && (Lineage != 0))
                {
                  Extinct(l, t);
                  alive--;
                }
              else
                {
                  // Transfer lineage variables to next time
                  for (chi = 0; chi < NCHAR; chi++)
                    {
                      Char_ltc[l][t+1][chi] = Char_ltc[l][t][chi];
                    }
                  Living_lt[l][t+1] = Living_lt[l][t];
                }
            }
        }

      t++;
      if ( ( (t == MAX_TIME ) && (strcmp(STOPPER, "time") == 0 ) ) || \
           ( (alive > LINEAGES ) && (strcmp(STOPPER, "lineages") == 0 ) ) || \
           (t  > BAILT ) ) stop = 1;
    }
  // need to reset MAX_TIME
  MAX_TIME = t;

  //fprintf(stderr, "time: %d, lineages: %d, ", t, alive);
  tout = t; aliveout = alive;

  // for those lineages still alive, kill em all! 
  //(thus giving them nodes & names)
  for (l = 0; l <= Lineage; l++)
    {
      Extant_l[l] = 0;
      if (Active_l[l] == 1)
        {
          Extinct(l, MAX_TIME);
          Extant_l[l] = 1;

        }
    }

  // New section: add current time extant traceback
  for (l = 0; l <= Lineage; l++)
    {
      if (Extant_l[l] == 1)
        {
          lroute = l;
          for (t = MAX_TIME; t >= 0; t--)
            {
              // Mark the censor route
              Censor_lt[lroute][t] = 1;
              if (Living_lt[lroute][t-1] == 0)
                {
                  lroute = LineageUp_l[lroute];
                }
            }
        }
    }


  // Now make lineage through time polots
  for (t = 0; t <= MAX_TIME; t++)
    {
      LttZ = 0;
      LttC = 0;
      for (chi = 0; chi < NCHAR; chi++)
        {
          CharMean[chi] = 0.0;
          CharSumXSqr[chi] = 0.0;
          CharCMean[chi] = 0.0;
          CharCSumXSqr[chi] = 0.0;
        }
      for (l = 0; l <= Lineage; l++)
        {
          if (Living_lt[l][t] == 1)
            {
              LttZ++;
              for (chi = 0; chi < NCHAR; chi++)
                {
                  CharMean[chi] += (float) Char_ltc[l][t][chi];
                  CharSumXSqr[chi] += (float) (Char_ltc[l][t][chi] * Char_ltc[l][t][chi]);
                }
            }
          if (Censor_lt[l][t] == 1)
            {
              LttC++;
              for (chi = 0; chi < NCHAR; chi++)
                {
                  CharCMean[chi] += (float) Char_ltc[l][t][chi];
                  CharCSumXSqr[chi] += (float) (Char_ltc[l][t][chi] * Char_ltc[l][t][chi]);
                }
            }
        }
      // Calculate Char Variance
      
      for (chi = 0; chi < NCHAR; chi++)
        {
          CharSumXSqr[chi] = ( CharSumXSqr[chi] - ( ( CharMean[chi] * CharMean[chi]) / \
                                                    (float) LttZ ) ) / (float) (LttZ - 1);
          CharMean[chi] = CharMean[chi] / (float) LttZ;
      
          CharCSumXSqr[chi] = ( CharCSumXSqr[chi] - ( ( CharCMean[chi] * CharCMean[chi]) / \
                                                      (float) LttC ) ) / (float) (LttC - 1);
          CharCMean[chi] = CharCMean[chi] / (float) LttC;
        }
      // output LTT
      if (OUT_MODE == 2) printf("%d\t%d\t%6.1f\t%6.1f\t%d\t%6.1f\t%6.1f\n", t, LttZ, \
                                CharMean[0], CharSumXSqr[0], LttC, CharCMean[0], CharCSumXSqr[0]);
            
      // Print out
    }

  if (OUT_MODE == 3) 
    {
      MakePhylo();
      fprintf(stderr, ", time: %d, lineages: %d\n", tout, aliveout);
      WriteTraits();
      DummySample();
    }
    
  return 1;
}

/* FUNCTION: Speciate() ---------------------------------------------------- */

void Speciate(int l, int t)
{
  int x, chi;
  // make a new node
  Node++;
  Name++;
  Time_n[Node] = t;
  Name_n[Node] = Node;
  Terminal_n[Node] = 0;

  // what's the node above it?
  NodeUp_n[Node] = NodeUp_l[l];
  
  // what's the BL?
  BlUp_n[Node] = t - Time_n[NodeUp_n[Node]];
  
  // The Char value?
  for (chi = 0; chi < NCHAR; chi++)
    {
      Char_nc[Node][chi] = Char_ltc[l][t][chi];
    }  
  // make new lineages (can add polytomies)
  for (x = 0; x <=1; x++)
    {
      Lineage++;
      Active_l[Lineage] = 1;
      Living_lt[Lineage][t+1] = 1;
      NodeUp_l[Lineage] = Node;
      LineageUp_l[Lineage] = l;
      for (chi = 0; chi < NCHAR; chi++)
        {
          Char_ltc[Lineage][t+1][chi] = Char_ltc[l][t][chi];
        }
    }
  
  // deactivate current lineage
  Active_l[l] = 0;
  Living_lt[l][t+1] = 0;
    
  // see it:
  Output(t);
}

/* FUNCTION: CharChange() ------------------------------------------------------ */

void CharChange(int l, int t)
{
  int change, i;
  // int amount[17] = {-3,-2,-1,-1,-1,0,0,0,0,0,0,0,1,1,1,2,3};
  // change = (int) ((17.0 * (float) random()) /    
  //          (float) (RAND_MAX+1.0));

  for(i = 0; i < NCHAR; i++)
    {
      change = (int) (((double) charch.n * (double) random()) / \
                      (double) (RAND_MAX+1.0));

      // printf("change = %d, n = %d\n", change, charch.n);
      //if (OUT_MODE == 0) printf("change = %d\n", change);
  
      if (TAPER)
        {
          Char_ltc[l][t][i] = Char_ltc[l][t][i] + \
            (charch.x[change] * exp( ( -1.0 * (float) t ) / TAPERFACT ) );
          // View it: printf("t %d   l %d   d  %f\n", t, l, (charch.x[change] * exp( ( -1.0 * (float) t ) / TAPERFACT ) ));
        }
      else
        {
          Char_ltc[l][t][i] = Char_ltc[l][t][i] + charch.x[change];
        }
      // record max and min
      if (CharMin[i] > Char_ltc[l][t][i]) CharMin[i] = Char_ltc[l][t][i];
      if (CharMax[i] < Char_ltc[l][t][i]) CharMax[i] = Char_ltc[l][t][i];       
    }
}

/* FUNCTION: Compete() ------------------------------------------------------ */

void Compete(int t)
{
  // just competes on trait 0 at moment
  // finds the minimum EcoDist to each taxon

  // vvv EDIT HERE FOR COMPETITION ALG vvv

  int x,y,z;
  for (y = 0; y <= Lineage; y++)
    {
      // reset the ecodist to account for recent extinctions
      EcoDist_l[y] = 10000;
      for (x = 0; x <= Lineage; x++)
        {
          if ((Active_l[x] == 1) && (Active_l[y] == 1) && (x != y))
            {
              z = (int) fabs((double)((Char_ltc[x][t][0] - Char_ltc[y][t][0])));
              if (EcoDist_l[y] > z) EcoDist_l[y] = z;
            }
        }
    }
  // ^^^                               ^^^
}

void Extinct(int l, int t)
{
  int chi;
  // make a new node
  Node++;   
  Time_n[Node] = t;
  Name++;
  Name_n[Node] = Node;
  for (chi = 0; chi < NCHAR; chi++)
    {
      Char_nc[Node][chi] = Char_ltc[l][t][chi];
    }
  if (t == MAX_TIME) Terminal_n[Node] = 1;
  else Terminal_n[Node] = 2;

  // what's the node above it?
  NodeUp_n[Node] = NodeUp_l[l];
  
  // what's the BL?
  BlUp_n[Node] = t - Time_n[NodeUp_n[Node]];

  // deactivate current lineage
  Active_l[l] = 0;

  // see it:
  Output(t);

}

void Output(int t)
{
  if (OUT_MODE == 0)
    printf("EVENT t %d\t\tn %d\tu %d\tb %d\t' %d\tc %f\n", t, Node, NodeUp_n[Node], \
           BlUp_n[Node], Name_n[Node], Char_nc[Node][0]);
  
  if (OUT_MODE == 1) printf("%d\t%d\t%d\t%d\t%f\t%d\n", Node, NodeUp_n[Node], \
                            BlUp_n[Node], Name_n[Node], Char_nc[Node][0], Time_n[Node]);

  
}

void MakePhylo()
{
  // needs to be final, so we can initialize the arrays based on Lineage results
  char tmp[10];
  int x;
  int *prunex;

  // node runs from 1 for root to Node

  prunex = ivector(0, Node-1);

  OutTree.nnodes = Node;
  OutTree.up = ivector(0, Node-1);
  OutTree.noat = ivector(0, Node-1);
  OutTree.depth = ivector(0, Node-1);
  OutTree.bl = vector(0, Node-1);
  OutTree.age = vector(0, Node-1);
  OutTree.taxon = cmatrix(0, Node-1, 0, 10);
  OutTree.ntaxa = Node;
  OutTree.depth[0] = 0;

  for (x = 1; x <= Node; x++)
    { 

      if (Terminal_n[x] == 1) prunex[x-1] = 1;
      else prunex[x-1] = 0;

      OutTree.up[x-1] = NodeUp_n[x]-1;
      OutTree.bl[x-1] = (float) BlUp_n[x];
      
      if (Terminal_n[x] == 0) sprintf(tmp, "node%d", Name_n[x]-1);
      if (Terminal_n[x] == 1) sprintf(tmp, "sp%d", Name_n[x]-1);
      if (Terminal_n[x] == 2) sprintf(tmp, "dead%d", Name_n[x]-1);

      strcpy(OutTree.taxon[x-1], tmp);

      OutTree.age[x-1] = (float) Time_n[x];
      if (x != 1)
        {
          OutTree.noat[NodeUp_n[x]-1]++;
          OutTree.depth[x-1] = OutTree.depth[NodeUp_n[x]-1]+1;
        }
    }

  //printf("nnodes = %d\n",OutTree.nnodes);

  //  for (x = 0; x < OutTree.nnodes; x++)
  //  { 
  //      printf("%d\t%d\t%d\t%f\t%d\t%s\n", x, OutTree.up[x], OutTree.depth[x], OutTree.bl[x], OutTree.noat[x], OutTree.taxon[x]);

  //  }

  if (PRUNED == 1 ) Fy2newRec(Prune(OutTree, prunex));
  else Fy2newRec(OutTree);

  // Balance
  if (PRUNED == 1 ) 
    {
      fprintf(stderr, "balance: %f", Balance(Prune(OutTree, prunex)));
    }
  else 
    {
      for (x = 0; x < OutTree.nnodes; x++) 
        {
          if (OutTree.noat[x]==0) prunex[x] = 1;
          else prunex[x] = 0;
        }
      fprintf(stderr, "balance: %f", Balance(Prune(OutTree, prunex)));
    }  

}

shift MakeChange(char brownian[11])
{
  //char brownian[11] = "7000400200";
  int x, y, z;
  int brownsum0 = 0;
  int brownsum1 = 0;
  char tmp[2];
  int brownvec[10];
  shift brown;

  // sum them up
  for (x = 0; x < 10; x++)
    {
      strncpy(tmp, &brownian[x],1);
      sscanf(tmp, "%d", &brownvec[x]); 
      if (x == 0) brownsum0 += brownvec[x];
      if (x > 0) brownsum1 += brownvec[x];
    }
  
  brown.n = brownsum1 + brownsum0 + brownsum1;

  //  printf("length = %d\n",  (brownsum1 + brownsum0 + brownsum1));

  brown.x = ivector(0, (brownsum1 + brownsum0 + brownsum1)-1);

  x = 0;
  for (y = 9; y>0; y--)
    {
      for (z = 0; z < brownvec[y];z++)
        {
          brown.x[x] = y * -1;
          x++;
        }
    }
  for (y = 0; y<10; y++)
    {
      for (z = 0; z < brownvec[y]; z++)
        {
          brown.x[x] = y ;
          x++;
        }
    }

  //  for (x = 0; x < brown.n; x++)
  //{
  //  printf("%d\t%d\n", x, brown.x[x]);
  //}

  return brown;

}

void WriteTraits()
{

  int n, chi;
  FILE *out_file;

  if ( (out_file = fopen("ecovolve.traits", "w") ) == NULL )
    {   
      fprintf(stderr, "Error: cannot write to file ecovolve.traits\n");
      exit(8);
    }

  fprintf(out_file, "type");
  for (chi = 0; chi < NCHAR; chi++)
    {
      fprintf(out_file, "\t3");
    }
  fprintf(out_file, "\n");
  fprintf(out_file, "name");
  for (chi = 0; chi < NCHAR; chi++)
    {
      fprintf(out_file, "\tecov%d", chi);
    }
  fprintf(out_file, "\n");

  for (n = 1; n <= Node; n++)
    {
      if (PRUNED == 1)
        {
          if (Terminal_n[n] == 1)
            {
              fprintf(out_file, "%s", OutTree.taxon[n-1]);
              for (chi = 0; chi < NCHAR; chi++)
                {
                  fprintf(out_file, "\t%f", Char_nc[n][chi]);
                }
              fprintf(out_file, "\n");
            }
        }
      else
        {
          if (Terminal_n[n] != 0)
            {
              fprintf(out_file, "%s", OutTree.taxon[n-1]);
              for (chi = 0; chi < NCHAR; chi++)
                {
                  fprintf(out_file, "\t%f", Char_nc[n][chi]);
                }
              fprintf(out_file, "\n");
            }
        }
    }
  fclose(out_file);

}

void DummySample()
{

  int n;
  FILE *out_file;

  if ( (out_file = fopen("ecovolve.sample", "w") ) == NULL )
    {   
      fprintf(stderr, "Error: cannot write to file ecovolve.sample\n");
      exit(8);
    }

  for (n = 0; n < Node; n++)
    {
      if (Terminal_n[n+1] == 1)
        {
          fprintf(out_file, "alive\t1\t%s\n", OutTree.taxon[n]);
        }
    }
  fclose(out_file);

}

float Balance(struct phylo A)
{
  int i, j, k, nterm;
  int *nsub;
  int dev;
  float sumBal = 0.0;
  
  nsub = ivector(0, A.nnodes);

  for (i = 0; i < A.nnodes; i++)
    {
      //      printf("noat %d  ldown %d   rsister %d\n", A.noat[i], A.ldown[i], A.rsister[i]) ;
      nsub[i] = 0;

      // for just the terminal
      if (A.noat[i] == 0)
        {
          j = i;
          while (j >= 0)
            {
              nsub[j]++;
              j = A.up[j];
            }
        }
    }
  
  k = 0; nterm=0;
  for (i = 0; i < A.nnodes; i++)
    {
      // just for internal nodes
      if(A.noat[i] != 0)
        {
          dev = abs( nsub[A.ldown[i]] - nsub[A.rsister[A.ldown[i]]] );
          // sumBal += (float) dev / (float) nsub[i]; // my metric
          sumBal += (float) dev;  // Colless
          k++;
          // printf("%d\tk:%d\tL:%d\tR:%d\tdev:%d\tsumBal:%f\n", i, k, nsub[A.ldown[i]],  nsub[A.rsister[A.ldown[i]]], dev, sumBal);
        }
      else nterm++;
    }

  free_ivector(nsub, 0, A.nnodes);

  // printf("Bal: %f\n", sumBal / (float) k);
  // printf("Bal: %f\n", sumBal); // Colless (non Heard)
  return (sumBal / ((((float) nterm - 1.0) * ((float) nterm - 2.0)) / 2)); // Heard 92


}
