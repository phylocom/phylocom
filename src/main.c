#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "phylocom.h"
#include "nrutil.h"

int main(int argc, char *argv[])
{
  // int i, argx, exactmean, exactruns;
  int i, argx;

  // dimension arrays; initialize globals
  UseAbund = 0;
  LowSig = 1;
  NoBL = 0;
  Debug = 0;
  Verbose = 0;
  RUNS = MAXRUNS;
  SWAPS = MAXSWAPS;
  BURNIN = 0;
  XVAR = 1;
  AOTOUT = 0;
  HILEVEL = MAXLEVEL;
  //default swap is "samples become random draws from phylogeny"
  SWAPMETHOD = 2;  
  strcpy(PhyloFile, "phylo");
  strcpy(SampleFile, "sample");
  strcpy(TraitFile, "traits");
  RNDPRUNEN = 5;
  RNDPRUNET = 3;
  MAKENODENAMES = 0;
  NULLTESTING = 0;
  Droptail = 0;
  FYOUT = 0;

  // initialise functions:

  // set the buffer to zero, so output is written directly to file:
  setbuf(stdout, NULL);
  //SWK note - changed the following from srandom to srand, otherwise
  //the use of rand is not initialized! Also set to seed from time
  //so that multiple runs of the program give different answers
  //for randomization tests
  srandom(time(NULL));  // do this once only


  // read the arguments
  if (argc > 1)
    {
      sscanf(argv[1], "%s", Method);
    }
  else
    {
      PrintWelcome();
      return 1;
    }

  // pick up the global switches:
  for (argx = 0; argx <  argc; argx++)
    {
      // runs
      if (strcmp(argv[argx], "-r") == 0)
	{
	  sscanf(argv[argx+1], "%d", &RUNS);
	  if (strcmp(Method, "rndprune") == 0) RNDPRUNET = RUNS;
	}

      // prunes
      if (strcmp(argv[argx],"-p") == 0)
	{
	  sscanf(argv[argx+1],"%d", &RNDPRUNEN);
	}

      // swaps
      if (strcmp(argv[argx],"-w") == 0)
	{
	  sscanf(argv[argx+1],"%d", &SWAPS);
	}

      // swaps
      if (strcmp(argv[argx],"-b") == 0)
	{
	  sscanf(argv[argx+1],"%ld", &BURNIN);
	}


      // swap method/model
      // defaults to 0 (shuffle phylogeny)
      if (strcmp(argv[argx],"-m") == 0)
	{
	  sscanf(argv[argx+1],"%d", &SWAPMETHOD);
	}

      // phylo file
      if (strcmp(argv[argx], "-f") == 0)
	{
	  sscanf(argv[argx+1], "%s", PhyloFile);
	}

      if (strcmp(argv[argx],"-n") == 0)
	{
	  MAKENODENAMES = 1;
	  NULLTESTING = 1;
	}

      if (strcmp(argv[argx],"-a") == 0)
	{
	  UseAbund = 1;
	}


      // sample file
      if (strcmp(argv[argx], "-s") == 0)
	{
	  sscanf(argv[argx+1], "%s", SampleFile);
	}

      // traits file
      if (strcmp(argv[argx], "-t") == 0)
	{
	  sscanf(argv[argx+1], "%s", TraitFile);
	}

      // Ignore branch lengths in node calculations and taxon distances (traits)
      //TODO implement this for other functions as well
      if (strcmp(argv[argx], "-e") == 0)
	{
	  NoBL = 1;
	}

      if (strcmp(argv[argx], "-y") == 0)
        {
          FYOUT = 1;
        }

      //Debug
      if (strcmp(argv[argx], "-d") == 0)
	{
	  Debug = 1;
	}

      //Verbose - switch to toggle verbose output
      if (strcmp(argv[argx], "-v") == 0)
	{
	  Verbose = 1;
	}

      // Select x variable for contrasts (default = 0)
      if (strcmp(argv[argx], "-x") == 0)
	{
	  sscanf(argv[argx+1], "%d", &XVAR);
	}
      // Subset aot output (default = 0 = all)
      if (strcmp(argv[argx], "-o") == 0)
	{
	  sscanf(argv[argx+1], "%d", &AOTOUT);
	}

  }


  // MAIN METHOD SWITCHING:

  // HELP
  if (strcmp(Method, "help") == 0)
    {
      PrintWelcome();
      return 1;
    }

  // Simple new2fy conversion
  else if (strcmp(Method, "new2fy") == 0)
    {
      phylo Intree = New2fy(PhyloFile);

      for (i = 0 ; i < Intree.nnodes; i++)
	{
	  printf("%d\t%d\t%d\t%d,%d,%d...\t%d\t%f\t%s\n",\
		 i, Intree.up[i], Intree.noat[i], Intree.down[i][0], \
		 Intree.down[i][1],Intree.down[i][2],Intree.depth[i], \
		 Intree.bl[i], Intree.taxon[i]);
	}

      // Fy2newRec(Intree);

      return 1;
    }

  else if (strcmp(Method, "makenex") == 0)
    {
      phylo IntreeV[0];
      IntreeV[0] = ReadPhylogeny(PhyloFile);
      WriteNexus( IntreeV , 1, \
			  ReadSample(SampleFile), 1 , \
			  ReadTraits(TraitFile), 1 );

      return 1;
    }

  else if (strcmp(Method, "naf") == 0)
    {
      MAKENODENAMES = 1;
      phylo IntreeV[0];
      IntreeV[0] = ReadPhylogeny(PhyloFile);
      NAF( IntreeV , ReadSample(SampleFile), ReadTraits(TraitFile));

      return 1;
    }

  else if (strcmp(Method, "new2nex") == 0)
    {
      //SWK Changed this to use new NewickToNexus function
      //SWK NewickToNexus doesn't require traits/sample file - just phylo
      //SWK Should add option to use New2Nex function which includes trait/sample info
      NewickToNexus(ReadPhylogeny(PhyloFile));
      return 1;
    }

  else if (strcmp(Method, "cleanphy" ) == 0)
    {
      CleanPhy( ReadPhylogeny(PhyloFile) );
      return 1;
    }

  else if (strcmp(Method, "pd") == 0)
    {
      PD( ReadPhylogeny(PhyloFile), ReadSample(SampleFile) );

      return 1;
    }

  else if (strcmp(Method, "comnode") == 0)
    {
      Comnode(ReadPhylogeny("tree1"),ReadPhylogeny("tree2"));
      return 1;
    }

  else if (strcmp(Method, "bladj") == 0)
    {
      Bladj(ReadPhylogeny(PhyloFile));
      return 1;
    }

  else if (strcmp(Method, "version") == 0)
    {
      printf("%s; %s\n", VERSION, SVNREV);
      return 1;
    }

  else if (strcmp(Method, "license") == 0)
    {
      License();
      return 1;
    }

  else if (strcmp(Method,"agenode") == 0)
    {
      phylo Intree = ReadPhylogeny(PhyloFile);
      Intree.age = vector(0, Intree.nnodes-1);
      AgeNodes(Intree);
      for (i = 0; i < Intree.nnodes; i++)
	{
	  if (Intree.noat[i] >0)
	    {
	      printf("%s\t%f\n", Intree.taxon[i], Intree.age[i]);
	    }
	}
      return 1;
    }

  if (strcmp(Method,"ageterm") == 0)
    {
      phylo Intree = New2fy(PhyloFile);
      Intree.age = vector(0, Intree.nnodes-1);
      AgeNodes(Intree);
      for (i = 0; i < Intree.nnodes; i++)
	{
	  if (Intree.noat[i] == 0)
	    {
	      printf("%s\t%f\n", Intree.taxon[i], Intree.age[Intree.up[i]]);
	    }
	}
      return 1;
    }

  else if (strcmp(Method, "phydist") == 0)
    {
      SimpleDist(ReadPhylogeny(PhyloFile));
      return 1;
    }

  else if (strcmp(Method, "phyvar") == 0)
    {
      PhyloVarCovar(ReadPhylogeny(PhyloFile));
      return 1;
    }

  else if (strcmp(Method, "comdist") == 0)
    {
      // ComDist(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      CommunityDistanceNull(ReadPhylogeny(PhyloFile), ReadSample(SampleFile), SWAPMETHOD, UseAbund, NULLTESTING);
      return 1;
    }

  else if ((strcmp(Method, "comdistnn") == 0) || \
	   (strcmp(Method, "comdistnt") == 0))
    {
      // ComDist(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      CommunityDistanceNNNull(ReadPhylogeny(PhyloFile), ReadSample(SampleFile), SWAPMETHOD, UseAbund, NULLTESTING);
      return 1;
    }

  else if (strcmp(Method, "comdistold") == 0)
    {
      // ComDist(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      ComDist(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

  else if (strcmp(Method, "comdistnnold") == 0)
    {
      ComDistNN(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

  else if (strcmp(Method, "icomdist") == 0)
    {
      IComDist(New2fy(PhyloFile), ReadSample(SampleFile));
      IComDistNN(New2fy(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

  else if (strcmp(Method, "vcomdist") == 0)
    {
      VComDist(New2fy(PhyloFile), ReadSample(SampleFile));
      VComDistNN(New2fy(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

  else if (strcmp(Method, "means") == 0)
    {
      VMeans(New2fy(PhyloFile));
      return 1;
    }

  else if (strcmp(Method, "nodesig") == 0)
    {
      TreeView = 0;
      NodeSig(New2fy(PhyloFile), ReadSample(SampleFile), 0, UseAbund);
      return 1;
    }
  else if (strcmp(Method, "nodesigl") == 0)
    {
      TreeView = 0;
      NodeSig(New2fy(PhyloFile), ReadSample(SampleFile), 1, UseAbund);
      return 1;
    }

  /* else if (strcmp(Method, "diversity") == 0)
    {
      PhyloDiversity(New2fy(PhyloFile), ReadSample(SampleFile));
      return 1;
    }
  */

  else if (strcmp(Method, "rao") == 0)
    {
      RaoDiversity(New2fy(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

  else if (strcmp(Method, "ltt") == 0)
    {
      Ltt( ReadPhylogeny(PhyloFile), ReadSample(SampleFile) );
      return 1;
    }

  else if (strcmp(Method, "lttr") == 0)
    {
      LttR(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

//  else if (strcmp(Method, "slide") == 0)
//    {
//      Slide();
//      return 1;
//    }

  // DDA:
//  else if (strcmp(Method, "newaot") == 0)
//    {
//      AOT(ReadTraits(TraitFile), ReadPhylogeny(PhyloFile), 0, XVAR);
//      return 1;
//    }

  else if (strcmp(Method, "aot") == 0)
    {
      PSigRun(ReadTraits(TraitFile), ReadPhylogeny(PhyloFile), 0, XVAR, AOTOUT);
      return 1;
    }

  else if (strcmp(Method, "aotf") == 0)
    {
      PSigRun(ReadTraits(TraitFile), ReadPhylogeny(PhyloFile), 1, XVAR, AOTOUT);
      return 1;
    }

  else if (strcmp(Method, "aotn") == 0)
    {
      PSigRun(ReadTraits(TraitFile), ReadPhylogeny(PhyloFile), 2, XVAR, AOTOUT);
      return 1;
    }
  else if (strcmp(Method, "rndprune") == 0)
    {
      RandPrune(ReadPhylogeny(PhyloFile), RNDPRUNEN, RNDPRUNET);
      return 1;
    }

  else if (strcmp(Method, "sampleprune") == 0)
    {
      SamplePrune(ReadPhylogeny(PhyloFile), ReadSample(SampleFile));
      return 1;
    }

  // SWK:
  else if (strcmp(Method, "comstruct") == 0)
    {
      ComStruct(ReadPhylogeny(PhyloFile),ReadSample(SampleFile),SWAPMETHOD,UseAbund);
      return 1;
    }

  else if (strcmp(Method, "swap") == 0)
    {
      OutputSwappedMatrix(ReadPhylogeny(PhyloFile),ReadSample(SampleFile),SWAPMETHOD);
      return 1;
    }

  else if (strcmp(Method, "comtrait") == 0)
    {
      ComTraitMetric(ReadSample(SampleFile),ReadTraits(TraitFile),SWAPMETHOD, XVAR);
      return 1;
    }

  // Incorrect input:
  else
    {
      printf("\n Oops! Command not recognized.  Try 'phylocom help' for suggestions.\n");
      return 1;
    }
}

