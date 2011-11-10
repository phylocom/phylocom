// traits.c;  

#include <stdio.h>  
#include <stdlib.h>  
#include <math.h>  
#include <string.h>  
#include "phylocom.h"  
#include "nrutil.h"  
#include "stats.h"

// ---------- PSIG RUN --------------------------------------------------------  

void PSigRun(traits T, phylo P, int Output, int xVar, int aotOut) {
  //PSigRun is a control loop for trait analyses; counters for significance   
  //under null models are initialized to zero, and then    
  //observed values are counted as one in first rep. For each loop (run  
  // in RUNS) RandArray is called to randomize trait values across tips, and  
  // then trait statistics are calculated.  When run=0, observed data are  
  // used without randomization. Print statements output results and  
  // significance for trait conservatism and independent contrasts  

  int i, j, k, cx, run, depth, node, dtx, dtxlow, dtxhi, swap, tmp, bin;
  int RunNum = 0; // set to 1 to track run # for big jobs
  //int Debug = 0; // set to 1 for verbose output during development
  int inclContrasts = 0;
  int printNullResults = 0;
  // if 1: Tip.mn for traits 2 and greater is calculated only if trait 1 = 1
  int subSmpl = 0;
  int showDOT = 0;
  int DOTuseSD = 1; //wtdNodeAge weighted by sdNodeChar vs. contrast value
  // if (NoBL == 0) showDOT = 0;
  int LowSig = 1; // 1 = one-tailed stats for low values, 0=1tailed for high
  int t = 0;
  int z = 1;
  div_t d;
  float tmpf1, tmpf2;

  int **RndArr;
  int *PTattach;

  int maxDepth;
  //phylo P; //internal short name to improve formatting 
  // BEGIN CAM EDIT
  char temp[MAXNOTELENGTH];
  struct phylo Out[T.ntraits];
  // END CAM EDIT

  //dtx = number of descendent taxa in a contrast
  // dtxlow = number in low set; dtxhi = num in high set
  float tmpf, charv;
  //  float ssBG, ssTot;
  float **sumTipChar;
  float **sumsqTipChar;
  float **meanTipChar;
  float **sdTipChar;
  float **sumsqTipCharG; //sum squares of trait deviations from grand mean

  float **sumNodeChar;
  float **sumsqNodeChar;
  float **nodeChar;
  float **meanNodeCharBin;

  float **sumWtdNChar;
  float **sumsqWtdNChar;
  float **wtdNChar;
  float **sumBLNChar;
  float **sdNodeChar;

  int **addedAtNode;
  int **addedAtNodeBin;

  float **ssTipChar;
  float **ssWithinTips;
  float **ssAmongTips;
  float **percVarDivTips;

  float **ssTipsVNode;
  float **ssAmongNodes;
  float **ssWithinNodes;
  float **percVarDivNodes;

  //float **sumsqAmong;
  //float **sumAmong;
  float **ssBG1;
  float *sumSDTree;
  float *sumsqSDTree;
  float *charMean;
  float *BLwts;
  float *Tbl; //transformed branch lengths for weighted node values
  int **binChar;
  int **binCharAssign;
  int **nodeStatus;
  int **nodeStatusCons;

  // vars for correlation between node depth and contrast/st dev
  float *depthCorr;
  float *meanSDTree;
  float *sd_SDTree;
  float *sumDepth;
  float *sumNodeSd;
  float *sumprDpNvr;
  float *sumsqDepth;
  float *sumsqNodeSd;
  float ss1, ss2, sp;
  int *internalNodes;

  // vars for independent contrasts
  float **dChar;
  float **dBL;
  int *ordChar1;
  float dCharLow[T.ntraits];
  float dWtLow;
  float dCharHi[T.ntraits];
  float dWtHi;
  float contrast[T.ntraits];
  float **Contrast;
  float **contrasts;
  float varC;
  float *VarC;
  float *c;
  int *varCRankLow;
  int *varCRankHigh;
  float **DCharHi;
  float **DCharLow;
  float sumCon[T.ntraits];
  float sumsqCon[T.ntraits];
  float sumprodCon[T.ntraits];
  float tmpsum, med;
  float pic_corr[T.ntraits*2];
  float signCont[T.ntraits*2];
  float signContCons[T.ntraits*2];
  float sumConCons[T.ntraits*2];
  float sumsqConCons[T.ntraits*2];
  float sdCon[T.ntraits*2];
  float sdConCons[T.ntraits*2];
  float *sdContrast;
  float *SDContrast;
  int ordTrt[T.ntraits]; // to specify which trait is x

  //for SvS and DOT tests  
  float **AbsContrast;
  float minAContrast[T.ntraits];
  float maxAContrast[T.ntraits];
  float midAContrast[T.ntraits];
  float wtdContAge[T.ntraits];
  float ageWt[T.ntraits];
  float aveAgeLargeContrasts[T.ntraits];
  int nLargeContrasts[T.ntraits];
  int q[T.ntraits];
  int quadrant[T.ntraits][3];
  float SvS[T.ntraits*2];
  float SvSObs[T.ntraits*2];//GLOBAL TO STORE OBSERVED VALUE
  float DOT[T.ntraits*2];
  float DOTObs[T.ntraits*2];//GLOBAL
  int SvSSig[T.ntraits*2];
  int DOTSig[T.ntraits*2];
  float wtdContAgeObs[T.ntraits];

  // trait conservatism
  float **MeanNodeChar;
  float **MeanTipChar;
  float **WtdNChar;
  float **SdNodeChar;
  float **SdTipChar;
  float **SumsqTipCharG;
  float **SSTipChar;
  float **SSTipsVNode;
  float **SSAmongNodes;
  float **SSWithinNodes;
  float **PercVarDiv;
  float DepthCorr[T.ntraits];
  float MeanSDTree[T.ntraits];
  float SD_SDTree[T.ntraits];
  float SignCont[T.ntraits];
  float SumCon[T.ntraits];
  int **TipsAtNode;
  int NSignCont = 0;

  // independent contrasts and correlations
  float PICcorr[T.ntraits*2];
  int nContrast = 0;
  int nSignCont, nSignContCons;

  //new globals
  int NSignContCons = 0;
  float SignContCons[T.ntraits*2];
  float SumConCons[T.ntraits*2];
  float SDCon[T.ntraits*2];
  float SDConCons[T.ntraits*2];

  // Significance counters
  int **MeanSigL;
  int **SdSigL;
  int **TipMeanSigL;
  int **TipSdSigL;
  int **MeanSigH;
  int **SdSigH;
  int **TipMeanSigH;
  int **TipSdSigH;
  int DepthCorrSigL[T.ntraits*2];
  int SDTreeSigL[T.ntraits*2];
  int SD_SDTreeSigL[T.ntraits*2];
  int SignContSigL[T.ntraits*2];
  int PICcorrSigL[T.ntraits*2];
  int SignContConsSigL[T.ntraits*2];
  int PICcorrSignConsL[T.ntraits*2];
  int DepthCorrSigH[T.ntraits*2];
  int SDTreeSigH[T.ntraits*2];
  int SD_SDTreeSigH[T.ntraits*2];
  int SignContSigH[T.ntraits*2];
  int PICcorrSigH[T.ntraits*2];
  int SignContConsSigH[T.ntraits*2];
  int PICcorrSignConsH[T.ntraits*2];

  // printf("PSig: finish declarations");

  // contrast x variable - default is 0 - add flag later for user
  // to specify another variable as x
  cx = 0;

  //Allocate vectors and matrices

  addedAtNode = imatrix(0, T.ntraits-1, 0, P.nnodes);
  nodeChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  wtdNChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumWtdNChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumsqWtdNChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumBLNChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sdNodeChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumNodeChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumsqNodeChar = matrix(0, T.ntraits-1, 0, P.nnodes);

  BLwts = vector(0, P.nnodes);
  Tbl = vector(0, P.nnodes);

  addedAtNodeBin = imatrix(0, T.ntraits-1, 0, P.nnodes);
  meanNodeCharBin = matrix(0, T.ntraits-1, 0, P.nnodes);
  binChar = imatrix(0, T.ntraits-1, 0, P.nnodes);
  binCharAssign = imatrix(0, T.ntraits-1, 0, P.nnodes);

  meanTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sdTipChar = matrix(0, T.ntraits, 0, P.nnodes);
  sumTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumsqTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumsqTipCharG = matrix(0, T.ntraits-1, 0, P.nnodes);
  ssTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  ssWithinTips = matrix(0, T.ntraits-1, 0, P.nnodes);
  ssAmongTips = matrix(0, T.ntraits-1, 0, P.nnodes);

  nodeStatus = imatrix(0, T.ntraits-1, 0, P.nnodes);
  nodeStatusCons = imatrix(0, T.ntraits-1, 0, P.nnodes);
  sumSDTree = vector(0, T.ntraits-1);
  ssBG1 = matrix(0, T.ntraits-1, 0, P.nnodes);
  percVarDivTips = matrix(0, T.ntraits-1, 0, P.nnodes);
  percVarDivNodes = matrix(0, T.ntraits-1, 0, P.nnodes);
  PercVarDiv = matrix(0, T.ntraits-1, 0, P.nnodes);
  sumsqSDTree = vector(0, T.ntraits-1);
  charMean = vector(0, T.ntraits-1);
  ssTipsVNode = matrix(0, T.ntraits-1, 0, P.nnodes);
  ssAmongNodes = matrix(0, T.ntraits-1, 0, P.nnodes);
  ssWithinNodes = matrix(0, T.ntraits-1, 0, P.nnodes);
  SSTipsVNode = matrix(0, T.ntraits-1, 0, P.nnodes);
  SSAmongNodes = matrix(0, T.ntraits-1, 0, P.nnodes);
  SSWithinNodes = matrix(0, T.ntraits-1, 0, P.nnodes);

  depthCorr = vector(0, T.ntraits-1);
  meanSDTree = vector(0, T.ntraits-1);
  sd_SDTree = vector(0, T.ntraits-1);
  sumDepth = vector(0, T.ntraits-1);
  sumNodeSd = vector(0, T.ntraits-1);
  sumprDpNvr = vector(0, T.ntraits-1);
  sumsqDepth = vector(0, T.ntraits-1);
  sumsqNodeSd = vector(0, T.ntraits-1);
  internalNodes = ivector(0, T.ntraits-1);

  sdContrast = vector(0, P.nnodes);
  SDContrast = vector(0, P.nnodes);
  Contrast = matrix(0, T.ntraits-1, 0, P.nnodes);
  contrasts = matrix(0, T.ntraits-1, 0, P.nnodes);  
  c = vector(0, P.nnodes - P.termtaxa - 1);
  VarC = vector(0, T.ntraits-1);
  varCRankHigh = ivector(0,P.nnodes);
  varCRankLow = ivector(0,P.nnodes);
  DCharHi = matrix(0, T.ntraits-1, 0, P.nnodes);
  DCharLow = matrix(0, T.ntraits-1, 0, P.nnodes);
  AbsContrast = matrix(0, T.ntraits-1, 0, P.nnodes);
  MeanNodeChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  MeanTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  WtdNChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  SumsqTipCharG = matrix(0, T.ntraits-1, 0, P.nnodes);
  SdNodeChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  SdTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  SSTipChar = matrix(0, T.ntraits-1, 0, P.nnodes);
  TipsAtNode = imatrix(0, T.ntraits-1, 0, P.nnodes);
  TipMeanSigL = imatrix(0, T.ntraits-1, 0, P.nnodes);
  TipSdSigL = imatrix(0, T.ntraits-1, 0, P.nnodes);
  MeanSigL = imatrix(0, T.ntraits-1, 0, P.nnodes);
  SdSigL = imatrix(0, T.ntraits-1, 0, P.nnodes);
  TipMeanSigH = imatrix(0, T.ntraits-1, 0, P.nnodes);
  TipSdSigH = imatrix(0, T.ntraits-1, 0, P.nnodes);
  MeanSigH = imatrix(0, T.ntraits-1, 0, P.nnodes);
  SdSigH = imatrix(0, T.ntraits-1, 0, P.nnodes);

  PTattach = ivector(0, P.nnodes);
  AttachPhyloToTraits(P, T, PTattach);

  // Finish allocating memory

  node = 0;// prevents warning about using node uninitialized    
  ordTrt[0] = xVar - 1;
  t = 0;
  for (i = 1; i <= T.ntraits-1; i++) {
    if (t == xVar-1)
      t++;
    ordTrt[i] = t;
    if (T.type[t] > 1)
      inclContrasts = 1;
    t++;
  }

  if (Debug) {
    for (k = 0; k <= T.ntraits-1; k++) {
      printf("%d %d\n", k, ordTrt[k]);
    }
  }
  // exit(EXIT_FAILURE);
  // printf("%d\n",P.nnodes);

  for (t = 0; t <= T.ntraits-1; t++) {
    for (node = 0; node < P.nnodes; node++) {
      MeanSigL[t][node] = 0;
      SdSigL[t][node] = 0;
      TipMeanSigL[t][node] = 0;
      TipSdSigL[t][node] = 0;
      MeanSigH[t][node] = 0;
      SdSigH[t][node] = 0;
      TipMeanSigH[t][node] = 0;
      TipSdSigH[t][node] = 0;
      BLwts[node] = 0.0;
      Tbl[node] = 0.0;
    }
    SDTreeSigL[t] = 0;
    SD_SDTreeSigL[t] = 0;
    DepthCorrSigL[t] = 0;
    SDTreeSigH[t] = 0;
    SD_SDTreeSigH[t] = 0;
    DepthCorrSigH[t] = 0;
    varCRankLow[t] = 0;
    varCRankHigh[t] = 0;    
    for (j = 0; j <= 1; j++) {
      PICcorrSignConsL[t+j*T.ntraits] = 0;
      PICcorrSigL[t+j*T.ntraits] = 0;
      SignContSigL[t+j*T.ntraits] = 0;
      SignContConsSigL[t+j*T.ntraits] = 0;
      PICcorrSignConsH[t+j*T.ntraits] = 0;
      PICcorrSigH[t+j*T.ntraits] = 0;
      SignContSigH[t+j*T.ntraits] = 0;
      SignContConsSigH[t+j*T.ntraits] = 0;
      SvSSig[t+j*T.ntraits] = 0;
      DOTSig[t+j*T.ntraits] = 0;
    }
  }

  // Calculate grand mean for each trait for use in
  // calculation of sums of squares by node, below
  for (k = 0; k < T.ntraits; k++) {
    t = ordTrt[k];
    charMean[t] = 0.0;
    for (z = 0; z < T.ntaxa; z++) {
      charMean[t] += T.tr[z][t];
    }
    charMean[t] = charMean[t] / T.ntaxa;
  }

  RndArr = imatrix(0, T.ntaxa, 0, T.ntraits-1);

  AgeNodes(P);

  //Assign maxdepth - could also do this at time of phylo loading?
  maxDepth=0;
  for (i = 0; i < P.nnodes; i++) {
    if (maxDepth <= P.depth[i])
      maxDepth = P.depth[i];
  }

  for (run = 0; run <= RUNS; run++) {
    d = div(run, 50);
    if ((RunNum) && (d.rem == 0)) {
      printf("run %d\n", run);
    }

    RandArrayT(RndArr, T.ntaxa, run, T.ntraits);
    //printf("Enter PSig");  
    // CUT HERE
    //      
    nContrast = 0;
    nSignCont = 0;
    nSignContCons = 0;

    for (t = 0; t <= T.ntraits-1; t++) {
      sumSDTree[t] = 0.0;
      sumsqSDTree[t] = 0.0;
      internalNodes[t] = 0;
      sumNodeSd[t] = 0.0;
      sumsqNodeSd[t] = 0.0;
      sumDepth[t] = 0.0;
      sumsqDepth[t] = 0.0;
      sumprDpNvr[t] = 0.0;
      sumCon[t] = 0.0;
      sumsqCon[t] = 0.0;
      sumprodCon[t] = 0.0;
      sumConCons[t] = 0.0;
      sumsqConCons[t] = 0.0;
      signCont[t] = 0.0;
      signContCons[t] = 0.0;
      minAContrast[t] = 0.0;
      maxAContrast[t] = 0.0;
      nLargeContrasts[t] = 0;
      aveAgeLargeContrasts[t] = 0.0;
      wtdContAge[t] = 0.0;
      ageWt[t] = 0.0;

      for (node = 0; node < P.nnodes; node++) {
        sumNodeChar[t][node] = 0.0;
        sumsqNodeChar[t][node] = 0.0;
        sumsqTipCharG[t][node] = 0.0;
        addedAtNode[t][node] = 0;
        addedAtNodeBin[t][node] = 0;
        wtdNChar[t][node] = 0.0;
        sumWtdNChar[t][node] = 0.0;
        sumsqWtdNChar[t][node] = 0.0;
        sumBLNChar[t][node] = 0.0;
        sdNodeChar[t][node] = 0.0;
        nodeChar[t][node] = 0.0;
        meanNodeCharBin[t][node] = 0.0;
        sumTipChar[t][node] = 0.0;
        sumsqTipChar[t][node] = 0.0;
        TipsAtNode[t][node] = 0;
        meanTipChar[t][node] = 0.0;
        sdTipChar[t][node] = 0.0;
        //sumsqAmong[t][node] = 0.0;
        //sumAmong[t][node] = 0.0;
        ssWithinNodes[t][node] = 0.0;
        ssAmongNodes[t][node] = 0.0;
        ssTipsVNode[t][node] = 0.0;
        ssBG1[t][node] = 0.0;
        ssTipChar[t][node] = 0.0;
        ssWithinTips[t][node] = 0.0;
        ssAmongTips[t][node] = 0.0;
        contrast[t] = 0.0;
        // Next four variables are flags for 
        // marking nodes for binary analysis
        binCharAssign[t][node] = 0;
        // 1: binary character has been assigned to this node
        binChar[t][node] = 0;
        // -1: mixed states above node, >= 0 are character states
        nodeStatus[t][node] = -1;
        // after algorithms below: 1 = node for contrast
        nodeStatusCons[t][node] = -1;
        // 0 = node below contrast, no more use of this branch
      }
    }

    for (depth = maxDepth; depth >= 0; depth--) {
      for (node = 0; node < P.nnodes; node++) {
        // Terminals: Pass up tree - distal nodes always hit before
        // any ancestral nodes above them
        if ((P.depth[node] == depth) && (P.noat[node] == 0)) {
          P.age[node] = 0;
          if (P.age[P.up[node]] == -1) {
            P.age[P.up[node]] = P.age[node] + P.bl[node];
          }
          for (k = 0; k <= T.ntraits-1; k++) {
            t = ordTrt[k];
            //Increment sums, sumsq, N for t vals at ancestral node
            //NodeChar vars pass node averages up
            //TipChar vars pass raw tip values up
            //WtdNChar - BL wtd node averages
            //charv: local var holding trait value
            charv = T.tr[ RndArr[ PTattach[node] ][ t ] ][ t ];
            // pass t values into nodeChar
            nodeChar[t][node] = charv;

            // summation of tip values
            if (subSmpl == 0) {
              sumTipChar[t][P.up[node]] += charv;
              sumsqTipChar[t][P.up[node]] += charv * charv;
              TipsAtNode[t][P.up[node]]++;
            } else {
              if (T.tr[PTattach[node]][0] == 1) {
                sumTipChar[t][P.up[node]] += charv;
                sumsqTipChar[t][P.up[node]] += charv * charv;
                TipsAtNode[t][P.up[node]]++;
              }
            }
            // summation of node values to up.node
            addedAtNode[t][P.up[node]]++;
            sumNodeChar[t][P.up[node]] += charv;
            sumsqNodeChar[t][P.up[node]] += charv * charv;
            sumWtdNChar[t][P.up[node]] += charv;
            sumsqWtdNChar[t][P.up[node]] += charv * charv;
            sumBLNChar[t][P.up[node]] += charv * (1/P.bl[node]);
            if ((run == 0) && (k == 0)) {
              BLwts[P.up[node]] += 1/P.bl[node];
              Tbl[node] = P.bl[node];
            }

            // Miscellany not being used right now
            //sumAmong[t][P.up[node]] += charv;
            //sumsqAmong[t][P.up[node]] += charv * charv;
            tmpf = (charv - charMean[t]);
            sumsqTipCharG[t][P.up[node]] += tmpf * tmpf;
            ssBG1[t][P.up[node]] += charv * charv;

            // If character is binary
            if (T.type[ordTrt[0]] == 0) {
              //Increment binary vars in meanNodeCharBin
              meanNodeCharBin[t][P.up[node]] += charv;
              addedAtNodeBin[t][P.up[node]]++;
              meanNodeCharBin[t][node] = charv;

              //Analysis of binary nodes for contrasts:
              if (inclContrasts) {
                if (k == 0) {
                  if ((0) && (run == 0)) {
                    printf(
                           "NSt: %d  %d  %d  %d  %d  %d\n",
                           node,
                           P.up[node],
                           binChar[ordTrt[0]][node],
                           binChar[ordTrt[0]][P.up[node]],
                           nodeStatus[ordTrt[0]][node],
                           nodeStatus[ordTrt[0]][P.up[node]]);
                  }

                  // Binary t distribution is not randomized
                  // This preserves the number of contrasts under
                  // the null model, so randomization doesn't mix
                  // up power vs. null pattern.
                  bin = (int) T.tr[PTattach[node]][ordTrt[0]];
                  binCharAssign[ordTrt[0]][node] = 1;
                  binChar[ordTrt[0]][node] = bin;

                  // If this is first descendant of up[node],
                  // pass binary state and set assign flag = 1
                  if (nodeStatus[ordTrt[0]][P.up[node]] < 1) {
                    if (binCharAssign[ordTrt[0]][P.up[node]]
                        == 0) {
                      binCharAssign[ordTrt[0]][P.up[node]]
                        = 1;
                      binChar[ordTrt[0]][P.up[node]]
                        = bin;
                      nodeStatus[ordTrt[0]][P.up[node]]
                        = -1;
                    }
                    // if not first descendant of up[node]
                    // then if binary state differs from previous
                    // binary passed to up[node], set contrast
                    // flag = 1 (nodeStatus)
                    // If first contrast moving up, then set 
                    // conservative contrast flag = 1 (node below
                    // terminal will always be okay for conservative
                    // contrast test.
                    else {
                      if (binChar[ordTrt[0]][node]
                          != binChar[ordTrt[0]][P.up[node]]) {
                        nodeStatus[ordTrt[0]][P.up[node]]
                          = 1;
                        if (nodeStatusCons[ordTrt[0]][P.up[node]]
                            != 0) {
                          nodeStatusCons[ordTrt[0]][P.up[node]]
                            = 1;
                        }
                      }
                    }
                  }

                  if ((0) && (run == 0)) {
                    printf(
                           "NSt: %d  %d  %d  %d  %d  %d\n",
                           node,
                           P.up[node],
                           binChar[ordTrt[0]][node],
                           binChar[ordTrt[0]][P.up[node]],
                           nodeStatus[ordTrt[0]][node],
                           nodeStatus[ordTrt[0]][P.up[node]]);
                  }

                }
              }
            }
          }
        }
        // END Terminals      
        // Interiors

        if ((P.depth[node] == depth) && (P.noat[node] != 0)) {
          if (P.age[P.up[node]] == -1) {
            P.age[P.up[node]] = P.age[node] + P.bl[node];
          }
          for (k = 0; k <= T.ntraits-1; k++) {
            t = ordTrt[k];
            // Calculate t mean, sd and SS for tip values
            meanTipChar[t][node] = sumTipChar[t][node]
              / (float) TipsAtNode[t][node];
            sdTipChar[t][node] = (sumsqTipChar[t][node]
                                  - ((sumTipChar[t][node] * sumTipChar[t][node])
                                     / (float) TipsAtNode[t][node]))
              / (float) (TipsAtNode[t][node] - 1);
            if (sdTipChar[t][node]<0)
              (sdTipChar[t][node]=0);
            ssTipChar[t][node] = sdTipChar[t][node]
              * (TipsAtNode[t][node] - 1);
            ssAmongTips[t][node] = ssTipChar[t][node]
              - ssWithinTips[t][node];
            percVarDivTips[t][node] = ssAmongTips[t][node]
              / ssTipChar[t][node];
            sdTipChar[t][node]
              = (float) (sqrt((double) sdTipChar[t][node]));

            // calculate node values
            if ((Debug) && (run == 0)) {
              printf("%d\t%f\t%f\t%f\t%d\n", node,
                     sumBLNChar[t][node], BLwts[node],
                     sumNodeChar[t][node], addedAtNode[t][node]);
            }
            if (NoBL == 0) {
              nodeChar[t][node] = sumBLNChar[t][node]
                / BLwts[node];
            } else {
              nodeChar[t][node] = sumNodeChar[t][node]
                / addedAtNode[t][node];
            }
            tmpf = nodeChar[t][node] - meanTipChar[t][node];
            ssTipsVNode[t][node] = ssTipChar[t][node]
              + TipsAtNode[t][node] * tmpf * tmpf;
            wtdNChar[t][node] = sumWtdNChar[t][node]
              / TipsAtNode[t][node];
            tmpf = wtdNChar[t][node] - nodeChar[t][node];
            ssAmongNodes[t][node]
              = (sumsqWtdNChar[t][node]
                 - ((sumWtdNChar[t][node]
                     * sumWtdNChar[t][node])
                    / TipsAtNode[t][node]));
            ssAmongNodes[t][node] = ssAmongNodes[t][node]
              + TipsAtNode[t][node] * tmpf * tmpf;
            percVarDivNodes[t][node] = ssAmongNodes[t][node]
              / (ssAmongNodes[t][node]
                 + ssWithinNodes[t][node]);
            sdNodeChar[t][node]
              = (sumsqNodeChar[t][node]
                 - ((sumNodeChar[t][node]
                     * sumNodeChar[t][node])
                    / addedAtNode[t][node]));
            if (sdNodeChar[t][node]<0) {
              (sdNodeChar[t][node]=0);
            }
            tmpf = nodeChar[t][node] - sumNodeChar[t][node]
              / addedAtNode[t][node];
            sdNodeChar[t][node] = sdNodeChar[t][node]
              + addedAtNode[t][node] * tmpf * tmpf;
            sdNodeChar[t][node] = sqrt(sdNodeChar[t][node]
                                       / (addedAtNode[t][node] - 1));

            if (DOTuseSD == 1) {
              wtdContAge[t] += sdNodeChar[t][node] * P.age[node];
              ageWt[t] += sdNodeChar[t][node];
            }

            if (printNullResults) {
              printf("%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n", run, t,
                     node, TipsAtNode[t][node],
                     meanTipChar[t][node], sdTipChar[t][node],
                     nodeChar[t][node], sdNodeChar[t][node]);
            }

            // Depth-contrast correlation calculations
            //sumDepth[t] += P.depth[node];
            //sumsqDepth[t] += P.depth[node] * P.depth[node];
            sumDepth[t] += P.age[node];
            sumsqDepth[t] += P.age[node] * P.age[node];
            sumNodeSd[t] += sdNodeChar[t][node];
            sumsqNodeSd[t] += sdNodeChar[t][node]
              * sdNodeChar[t][node];
            //sumprDpNvr[t] += P.depth[node] * sdNodeChar[t][node];
            sumprDpNvr[t] += P.age[node] * sdNodeChar[t][node];

            // if not the root, increment node and tip variables
            // for trait.
            if (node != 0) {
              addedAtNode[t][P.up[node]]++;
              if ((run == 0) && (k == 0)) {
                Tbl[node] = P.bl[node] + 1/BLwts[node];
              }
              sumBLNChar[t][P.up[node]] += nodeChar[t][node]
                / Tbl[node];
              if ((run == 0) && (k == 0)) {
                BLwts[P.up[node]] += 1/Tbl[node];
              }
              sumNodeChar[t][P.up[node]] += nodeChar[t][node];
              sumsqNodeChar[t][P.up[node]] += (nodeChar[t][node]
                                               * nodeChar[t][node]);
              sumWtdNChar[t][P.up[node]] += nodeChar[t][node]
                * TipsAtNode[t][node];
              sumsqWtdNChar[t][P.up[node]] += (nodeChar[t][node]
                                               * nodeChar[t][node] * TipsAtNode[t][node]);
              ssWithinNodes[t][P.up[node]]
                += ssTipsVNode[t][node];

              sumTipChar[t][P.up[node]] += sumTipChar[t][node];
              sumsqTipChar[t][P.up[node]]
                += sumsqTipChar[t][node];
              TipsAtNode[t][P.up[node]] += TipsAtNode[t][node];
              ssWithinTips[t][P.up[node]] += ssTipChar[t][node];
              sumsqTipCharG[t][P.up[node]]
                += sumsqTipCharG[t][node];
              ssBG1[t][P.up[node]] += (sumTipChar[t][node]
                                       * sumTipChar[t][node])
                / TipsAtNode[t][node];
            }

            if (T.type[ordTrt[0]] > 1) {
              nodeStatus[ordTrt[0]][node] = 1;
            }

            // If trait 0 is binary              
            if (inclContrasts) {
              if (T.type[ordTrt[0]] == 0) {
                //meanNodeCharBin will be binary if all values passed
                //up were the same.
                meanNodeCharBin[t][node]
                  = meanNodeCharBin[t][node]
                  /addedAtNodeBin[t][node];
                //If all binary vals passed up are the same
                if ((nodeStatus[ordTrt[0]][node] == -1)
                    && (node != 0)) // check for root 
                  {
                    meanNodeCharBin[t][P.up[node]]
                      += meanNodeCharBin[t][node];
                    addedAtNodeBin[t][P.up[node]]++;
                  }
                //For binary char, analyze node status
                if (k == 0) {
                  if ((Debug) && (run == 0)) {
                    printf(
                           "NSi: %d  %d  %d  %d  %d  %d\n",
                           node,
                           P.up[node],
                           binChar[ordTrt[0]][node],
                           binChar[ordTrt[0]][P.up[node]],
                           nodeStatus[ordTrt[0]][node],
                           nodeStatus[ordTrt[0]][P.up[node]]);
                  }
                  //If this is a conservative contrast node, or if
                  //any node below was one, then set flag = 0, so this
                  //node will not be used for conservative contrasts.
                  if (nodeStatusCons[ordTrt[0]][node] >= 0) {
                    nodeStatusCons[ordTrt[0]][P.up[node]]
                      = 0;
                  }
                  //If all binaries the same so far, then pass
                  //binary values up
                  if (nodeStatus[ordTrt[0]][node] == -1) {
                    if (binCharAssign[ordTrt[0]][P.up[node]]
                        == 0) {
                      binCharAssign[ordTrt[0]][P.up[node]]
                        = 1;
                      binChar[ordTrt[0]][P.up[node]]
                        = binChar[ordTrt[0]][node];
                      nodeStatus[ordTrt[0]][P.up[node]]
                        = -1;
                    }
                    //if there are contrasting binaries at this level
                    //then set contrast flag = 1; if this is first
                    //contrast on this path, then set conservative
                    //contrast flag = 1.
                    else {
                      if (binChar[ordTrt[0]][P.up[node]]
                          != binChar[ordTrt[0]][node]) {
                        nodeStatus[ordTrt[0]][P.up[node]]
                          = 1;
                        if (nodeStatusCons[ordTrt[0]][P.up[node]]
                            != 0) {
                          nodeStatusCons[ordTrt[0]][P.up[node]]
                            = 1;
                        }
                      }
                    }
                  }
                  //If contrast has already been calculated at this
                  //level or above, and if this is first branch
                  //to look up, then set contrast flag = 0
                  else {
                    if ((nodeStatus[ordTrt[0]][P.up[node]]
                         == -1)
                        && (binCharAssign[ordTrt[0]][P.up[node]]
                            == 0)) {
                      nodeStatus[ordTrt[0]][P.up[node]]
                        = 0;
                    }
                  }
                }

                if ((Debug) && (run == 0)) {
                  printf("NSi: %d  %d  %d  %d  %d  %d\n",
                         node, P.up[node],
                         binChar[ordTrt[0]][node],
                         binChar[ordTrt[0]][P.up[node]],
                         nodeStatus[ordTrt[0]][node],
                         nodeStatus[ordTrt[0]][P.up[node]]);
                }
              }
            }
            //End binary algorithms

            //Increment treewise sum of sd of divergences by node
            sumSDTree[t] += sdNodeChar[t][node];
            sumsqSDTree[t] += sdNodeChar[t][node]
              * sdNodeChar[t][node];
            internalNodes[t]++;

            // Pass to Globals for observed
            if (run == 0) {
              WtdNChar[t][node] = wtdNChar[t][node];
              MeanNodeChar[t][node] = nodeChar[t][node];
              SdNodeChar[t][node] = sdNodeChar[t][node];
              MeanTipChar[t][node] = meanTipChar[t][node];
              SdTipChar[t][node] = sdTipChar[t][node];
              SumsqTipCharG[t][node] = sumsqTipCharG[t][node];
              PercVarDiv[t][node] = percVarDivTips[t][node];
              SSTipChar[t][node] = ssTipChar[t][node];
              SSTipsVNode[t][node] = ssTipsVNode[t][node];
              SSAmongNodes[t][node] = ssAmongNodes[t][node];
              SSWithinNodes[t][node] = ssWithinNodes[t][node];
            }
            //Compare current values with observed (in globals)
            //and increment significance counters. Comparison is
            //made in first run to count observed as 1.

            if (nodeChar[t][node] <= MeanNodeChar[t][node]) {
              MeanSigL[t][node]++;
            }
            if (sdNodeChar[t][node] <= SdNodeChar[t][node]) {
              SdSigL[t][node]++;
            }
            if (meanTipChar[t][node] <= MeanTipChar[t][node]) {
              TipMeanSigL[t][node]++;
            }
            if (sdTipChar[t][node] <= SdTipChar[t][node]) {
              TipSdSigL[t][node]++;
            }
            if (nodeChar[t][node] >= MeanNodeChar[t][node]) {
              MeanSigH[t][node]++;
            }
            if (sdNodeChar[t][node] >= SdNodeChar[t][node]) {
              SdSigH[t][node]++;
            }
            if (meanTipChar[t][node] >= MeanTipChar[t][node]) {
              TipMeanSigH[t][node]++;
            }
            if (sdTipChar[t][node] >= SdTipChar[t][node]) {
              TipSdSigH[t][node]++;
            }
          }
        }
      }
    }

    for (k = 0; k <= T.ntraits-1; k++) {
      t = ordTrt[k];
      //Significance at root set to RUNS+1 - should be ignored in results
      TipMeanSigL[t][0] = RUNS+1;
      TipSdSigL[t][0] = RUNS+1;
      TipMeanSigH[t][0] = RUNS+1;
      TipSdSigH[t][0] = RUNS+1;
      //Calculation of treewise mean and sd of divergence sds
      meanSDTree[t] = sumSDTree[t] / (float) internalNodes[t];
      sd_SDTree[t]
        = sqrt((sumsqSDTree[t] - ((sumSDTree[t] * sumSDTree[t])
                                  / (float) internalNodes[t]))
               / (float) (internalNodes[t] - 1));
      //Depth correlation calculation 
      ss1 = sumsqDepth[t] - (sumDepth[t] * sumDepth[t])
        / internalNodes[t];
      ss2 = sumsqNodeSd[t] - (sumNodeSd[t] * sumNodeSd[t])
        / internalNodes[t];
      sp = sumprDpNvr[t] - (sumDepth[t] * sumNodeSd[t])
        / internalNodes[t];
      depthCorr[t] = sp / sqrt(ss1*ss2);

      //Pass to globals          
      if (run == 0) {
        MeanSDTree[t] = meanSDTree[t];
        SD_SDTree[t] = sd_SDTree[t];
        DepthCorr[t] = depthCorr[t];
      }

      // compare whole tree measures to 
      if (meanSDTree[t] <= MeanSDTree[t]) {
        SDTreeSigL[t]++;
      }
      if (sd_SDTree[t] <= SD_SDTree[t]) {
        SD_SDTreeSigL[t]++;
      }
      if (depthCorr[t] <= DepthCorr[t]) {
        DepthCorrSigL[t]++;
      }
      if (meanSDTree[t] >= MeanSDTree[t]) {
        SDTreeSigH[t]++;
      }
      if (sd_SDTree[t] >= SD_SDTree[t]) {
        SD_SDTreeSigH[t]++;
      }
      if (depthCorr[t] >= DepthCorr[t]) {
        DepthCorrSigH[t]++;
      }
      //see results from each randomization
      if (0)
        if (k == 0) {
          printf("%d\t%f\t%f\t%f\t%f\n", run, meanTipChar[t][0],
                 meanSDTree[t], sd_SDTree[t], depthCorr[t]);
        }
    }

    if ((Debug) && (run == 0)) {
      for (node = 0; node < P.nnodes; node++) {
        printf("%d\t%d\t%d\t%d\t%f\n", node, binChar[ordTrt[0]][node],
               nodeStatus[ordTrt[0]][node],
               nodeStatusCons[ordTrt[0]][node],
               meanNodeCharBin[1][node]);
      }
    }

    // INDEPENDENT CONTRASTS


    // Pagel (1992 J Theor Biol) contrast algorithm: 
    // trait values for both characters are ordered, and then
    // divided on median value for char 1 to create 1 contrast 
    if (inclContrasts) {
      for (node = 0; node < P.nnodes; node++) {
        if (P.noat[node] != 0)//If not terminal node
          {
            if ((Debug) && (run == 0))
              printf("entering contrast loop %d\t%d\n", node,
                     P.noat[node]);
            ordChar1 = ivector(0, P.noat[node]-1);
            dChar = matrix(0, T.ntraits-1, 0, P.noat[node]-1);
            dBL = matrix(0, T.ntraits-1, 0, P.noat[node]-1);
            for (i = 0; i <= P.noat[node] - 1; i++) // Initialize
              {
                ordChar1[i] = 0;//sort array
                for (t = 0; t <= T.ntraits - 1; t++) {
                  dChar[t][i] = 0.0;
                }
              }

            // Pass trait values for Down nodes into matrix for this node
            // HERE: may need to add binary dChar vector...
            dtx = 0;
            for (i = 1; i < P.nnodes; i++) {
              //if (0) 
              //  {
              //    printf("%d\t%d\n",i,P.up[i]);
              //  }
              if (P.up[i] == node) {
                //if (0) 
                //  {
                //    printf("MATCH");
                //  }
                for (k = 0; k <= T.ntraits-1; k++) {
                  t = ordTrt[k];
                  if (T.type[ordTrt[0]] == 0) {
                    //only pass binary value up if down node vals are
                    //all the same
                    if (nodeStatus[ordTrt[0]][i] == -1) {
                      if (k == 0) {
                        dtx++;
                        dChar[t][dtx-1]
                          = (float) binChar[t][i];
                      }
                      //passing meanNodeCharBin treats 
                      //chars 1 and up as normal characters.
                      else {
                        dChar[t][dtx-1]
                          = meanNodeCharBin[t][i];
                      }
                    }
                  } else {
                    if (k == 0)
                      dtx++;
                    dChar[t][dtx-1] = nodeChar[t][i];
                    dBL[t][dtx-1] = Tbl[i];
                  }
                }
                ordChar1[dtx-1] = dtx-1;
              }
            }

            // rank values for character 1 (t = 0); bubble sort
            do {
              swap = 0;
              for (i = 0; i < dtx-1; i++) {
                if (dChar[ordTrt[0]][ordChar1[i]]
                    > dChar[ordTrt[0]][ordChar1[i+1]]) {
                  tmp = ordChar1[i];
                  ordChar1[i] = ordChar1[i+1];
                  ordChar1[i+1] = tmp;
                  swap = 1;
                }
              }
            } while (swap == 1);

            if ((Debug) && (run == 0)) {
              for (i = 0; i < dtx; i++) {
                printf("dtx sort: %d\t%f\t%f\t%f\t%f\n", node,
                       dChar[ordTrt[0]][ordChar1[i]],
                       dChar[ordTrt[1]][ordChar1[i]],
                       dBL[ordTrt[0]][ordChar1[i]],
                       dBL[ordTrt[1]][ordChar1[i]]);
              }
            }

            //If character is binary, then split has to be between the two sets
            //of values, which won't be at the midpoint of dtx; this section
            //calculates where split is.
            if (T.type[ordTrt[0]] == 0) {
              bin = (int) dChar[ordTrt[0]][ordChar1[0]];
              dtxlow = 0;
              dtxhi = 0;
              for (i = 0; i < dtx; i++) {
                if ((int) dChar[ordTrt[0]][ordChar1[i]] == bin) {
                  dtxlow++;
                } else {
                  dtxhi++;
                }
              }
              if ((run == 0) && (dtx != dtxlow + dtxhi))
                printf("dtx sum error");
            } else //Here for continuous characters, splitting at midpoint
              {
                if (dtx%2 == 0) {
                  dtxlow = dtx/2;
                  dtxhi = dtx/2;
                } else {
                  // odd number of traits in polytomy, 
                  // decide where to put median value
                  tmpsum = 0.0;
                  for (i = 0; i < dtx; i++) {
                    tmpsum += dChar[ordTrt[0]][ordChar1[i]];
                  }
                  tmpsum = tmpsum/dtx;
                  med = dChar[ordTrt[0]][ordChar1[(int) ((float) (dtx-1)/2.0)]];
                  if (med < tmpsum) {
                    dtxlow = (dtx+1)/2;
                    dtxhi = (dtx-1)/2;
                  } else {
                    dtxlow = (dtx-1)/2;
                    dtxhi = (dtx+1)/2;
                  }
                }
              }
            if ((Debug) && (run == 0)) {
              printf("%d\t%d\t%d\n", dtx, dtxlow, dtxhi);
            }

            //Sum trait values in low and high halves and calculate means
            for (k = 0; k <= T.ntraits-1; k++) {
              t = ordTrt[k];
              {
                dCharLow[t] = 0.0;
                dCharHi[t] = 0.0;
              }
              dWtLow = 0.0;
              dWtHi = 0.0;
              if (T.type[ordTrt[0]] == 0) {
                for (i = 0; i < dtxlow; i++) {
                  dCharLow[t] += dChar[t][ordChar1[i]];
                }
                for (i = dtxlow; i <= dtx-1; i++) {
                  dCharHi[t] += dChar[t][ordChar1[i]];
                }
                dCharLow[t] = dCharLow[t]/dtxlow;
                dCharHi[t] = dCharHi[t]/dtxhi;
                contrast[t] = dCharHi[t] - dCharLow[t];
              } else {
                if (NoBL == 0) {
                  for (i = 0; i < dtxlow; i++) {
                    dCharLow[t] += dChar[t][ordChar1[i]]
                      / dBL[t][ordChar1[i]];
                    dWtLow += 1 / dBL[t][ordChar1[i]];
                  }
                  for (i = dtxlow; i <= dtx-1; i++) {
                    dCharHi[t] += dChar[t][ordChar1[i]]
                      / dBL[t][ordChar1[i]];
                    dWtHi += 1 / dBL[t][ordChar1[i]];
                  }
                } else {
                  for (i = 0; i < dtxlow; i++) {
                    dCharLow[t] += dChar[t][ordChar1[i]];
                    dWtLow++;
                  }
                  for (i = dtxlow; i <= dtx-1; i++) {
                    dCharHi[t] += dChar[t][ordChar1[i]];
                    dWtHi++;
                  }
                }
                // check harmonic means for BL!!! - is it N/summation?
                dCharLow[t] = dCharLow[t]/dWtLow;
                dCharHi[t] = dCharHi[t]/dWtHi;
                dWtLow = dtxlow / dWtLow;
                dWtHi = dtxhi / dWtHi;
                if (NoBL == 0) {
                  sdContrast[node] = sqrt(dWtLow + dWtHi);
                  contrast[t] = (dCharHi[t] - dCharLow[t])
                    / sdContrast[node];
                } else {
                  sdContrast[node] = 0.0;
                  contrast[t] = (dCharHi[t] - dCharLow[t]);
                  if (DOTuseSD == 0) {
                    wtdContAge[t] += fabs(contrast[t])
                      * P.age[node];
                    ageWt[t] += fabs(contrast[t]);
                  }
                }
                if ((Debug) && (run == 0)) {
                  printf("CONT:\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n",
                         node, t, dCharLow[t], dCharHi[t],
                         contrast[t], dWtLow, dWtHi);
                }
              }

              if (run == 0) {
                DCharHi[t][node] = dCharHi[t];
                DCharLow[t][node] = dCharLow[t];
                Contrast[t][node] = contrast[t];
                if (Debug) printf("CONTRAST\t%f\n",contrast[t]);
                SDContrast[node] = sdContrast[node];
                AbsContrast[t][node] = fabs(contrast[t]);
                if (minAContrast[t] == 0)
                  minAContrast[t] = fabs(contrast[t]);
                if (fabs(contrast[t]) < minAContrast[t]) {
                  minAContrast[t] = fabs(contrast[t]);
                }
                if (fabs(contrast[t]) > maxAContrast[t]) {
                  maxAContrast[t] = fabs(contrast[t]);
                }
              }
                        
              //for signal calculations
              contrasts[t][node] = contrast[t];

              //For binary trait, if contrast node increment signCont by
              //1 if contrast positive, and by 1/2 if contrast = 0
              if (T.type[ordTrt[0]] == 0) // trait 1 is binary
                {
                  if (nodeStatus[ordTrt[0]][node] == 1) {
                    //if ((run == 0)&&(t == 1)) printf("enter binary contrast calc\n");
                    if (k == 1)
                      nSignCont++;
                    if (contrast[t] > 0)
                      signCont[t] += 1;
                    if (contrast[t] == 0)
                      signCont[t] += 0.5;
                    sumCon[t] += contrast[t];
                    sumsqCon[t] += contrast[t] * contrast [t];
                  }
                  //Increment sums for conservative contrasts
                  if (nodeStatusCons[ordTrt[0]][node] == 1) {
                    //if ((run == 0)&&(t == 1)) printf("enter binary contrast calc\n");
                    if (k == 1)
                      nSignContCons++;
                    if (contrast[t] > 0)
                      signContCons[t] += 1;
                    if (contrast[t] == 0)
                      signContCons[t] += 0.5;
                    sumConCons[t] += contrast[t];
                    sumsqConCons[t] += contrast[t] * contrast[t];
                  }

                }
              // trait 1 is ordered multistate or continuous
              else if (T.type[ordTrt[0]] > 1) {
                sumCon[t] += abs(contrast[t]);
                sumsqCon[t] += contrast[t] * contrast [t];
                if (k != 0) {
                  if (k == 1)
                    nContrast++;
                  sumprodCon[t] += contrast[ordTrt[0]]
                    * contrast[t];
                  if (k == 1)
                    nSignCont++;
                  if (contrast[t] > 0)
                    signCont[t] += 1;
                  if (contrast[t] == 0)
                    signCont[t] += 0.5;
                }
              }
            }

            free_ivector(ordChar1, 0, P.noat[node]);
            free_matrix(dChar, 0, T.ntraits-1, 0, P.noat[node]);
            free_matrix(dBL, 0, T.ntraits-1, 0, P.noat[node]);

            if ((Debug) && (run == 0)) {
              printf("%d\t%d\t%f\t%f\n", nodeStatus[ordTrt[0]][node],
                     nodeStatusCons[ordTrt[0]][node], dCharLow[0],
                     dCharLow[1]);
              printf("%d\t%d\t%f\t%f\n", nodeStatus[ordTrt[0]][node],
                     nodeStatusCons[ordTrt[0]][node], dCharHi[0],
                     dCharHi[1]);
              printf("%d\t%d\t%f\t%f\n", nodeStatus[ordTrt[0]][node],
                     nodeStatusCons[ordTrt[0]][node], contrast[0],
                     contrast[1]);
            }

            if ((Debug) && (run == 0)) {
              printf("\n");
            }
          }
      }

      //SvS and DOT tests
      for (k = 0; k <= T.ntraits - 1; k++) {
        t = ordTrt[k];
        quadrant[t][0] = 0;
        quadrant[t][1] = 0;
        quadrant[t][2] = 0;
        midAContrast[t] = (maxAContrast[t] + minAContrast[t]) / 2;
        if (run == -1)
          printf("%f\t%f\n", minAContrast[0], minAContrast[1]);
        wtdContAge[t] = wtdContAge[t] / ageWt[t];
        if (run == 0)
          wtdContAgeObs[t] = wtdContAge[t];
      }
      //For node loop here - calculating quadrant scores and mean ages
      for (node = 0; node < P.nnodes; node++) {
        if (P.noat[node] != 0) // Internal nodes only
          {
            for (k = 0; k <= T.ntraits - 1; k++) {
              t = ordTrt[k];
              q[t] = 0;
              if (AbsContrast[t][node] >= midAContrast[t]) {
                q[t] = 1;
                nLargeContrasts[t]++;
                aveAgeLargeContrasts[t] += P.age[node];
                if (run == -1) {
                  printf("t,node,q,age %d\t%d\t%f\n", t, node,
                         P.age[node]);
                }
              }
              if (k > 0) {
                quadrant[t][q[0]+q[t]]++;
              }
            }
          }
      }

      //Calculate SvS and DOT
      aveAgeLargeContrasts[ordTrt[0]] = aveAgeLargeContrasts[ordTrt[0]]
        / nLargeContrasts[ordTrt[0]];
      for (k = 0; k <= T.ntraits - 1; k++) {
        t = ordTrt[k];
        SvS[t] = (float) quadrant[t][1]/ (float) (quadrant[t][1]
                                                  + quadrant[t][2]);
        aveAgeLargeContrasts[t] = aveAgeLargeContrasts[t]
          / nLargeContrasts[t];
        DOT[t] = wtdContAge[ordTrt[0]] - wtdContAge[t];
        if (run == 0) {
          SvSObs[t] = SvS[t];
          DOTObs[t] = DOT[t];
        }
        if (LowSig) {
          if (SvS[t] <= SvSObs[t])
            SvSSig[t]++;
          if (DOT[t] <= DOTObs[t])
            DOTSig[t]++;
        } else {
          if (SvS[t] >= SvSObs[t])
            SvSSig[t]++;
          if (DOT[t] >= DOTObs[t])
            DOTSig[t]++;
        }
      }
      if (run == -1) {
        printf("quadrants 1, 2 %d\t%d\n", quadrant[1][1],
               quadrant[1][2]);
        printf("Age1 %f\n", aveAgeLargeContrasts[ordTrt[0]]);
        printf("Age1 %f\n", aveAgeLargeContrasts[1]);
        printf("SvS %f, DOT %f\n", SvS[1], DOT[1]);
      }
      //End SvS and DOT tests        


      //Pass sign test to globals for comparison with nulls
      for (k = 1; k <= T.ntraits - 1; k++) {
        t = ordTrt[k];
        if (run == 0) {
          SignCont[t] = signCont[t];
          SignContCons[t] = signContCons[t];
          NSignCont = nSignCont;
          NSignContCons = nSignContCons;
        }
        if (signCont[t] <= SignCont[t]) {
          SignContSigL[t]++;
        }
        if (signContCons[t] <= SignContCons[t]) {
          SignContConsSigL[t]++;
        }
        if (signCont[t] >= SignCont[t]) {
          SignContSigH[t]++;
        }
        if (signContCons[t] >= SignContCons[t]) {
          SignContConsSigH[t]++;
        }
      }

      //Pass contrast correlations to globals, and compare with nulls to
      //increment significance counters
      if (T.type[ordTrt[0]] == 0) {
        for (k = 1; k <= T.ntraits-1; k++) {
          t = ordTrt[k];
          sdCon[t] = sumsqCon[t] - ((sumCon[t] * sumCon[t])
                                    / nSignCont);
          sdCon[t] = sqrt(sdCon[t] / (nSignCont - 1));
          sdConCons[t] = sumsqConCons[t] - ((sumConCons[t]
                                             * sumConCons[t]) / nSignContCons);
          sdConCons[t] = sqrt(sdConCons[t] / (nSignContCons - 1));
          sumCon[t] = sumCon[t]/nSignCont;
          sumConCons[t] = sumConCons[t]/nSignContCons;
          //printf("%f\n",sumCon[t]);
          if (run == 0) {
            SumCon[t] = sumCon[t];
            SumConCons[t] = sumConCons[t];
            SDCon[t] = sdCon[t];
            SDConCons[t] = sdConCons[t];
          }
          if (sumCon[t] <= SumCon[t]) {
            PICcorrSigL[t]++;
          }
          if (sumConCons[t] <= SumConCons[t]) {
            PICcorrSignConsL[t]++;
          }
          if (sumCon[t] >= SumCon[t]) {
            PICcorrSigH[t]++;
          }
          if (sumConCons[t] >= SumConCons[t]) {
            PICcorrSignConsH[t]++;
          }
          //printf("%f\n", PICcorr[t]);
        }
      } else {
        for (k = 1; k <= T.ntraits-1; k++) {
          t = ordTrt[k];
          pic_corr[t] = sumprodCon[t] / sqrt(sumsqCon[ordTrt[0]]
                                             * sumsqCon[t]);
          if (run == 0) {
            PICcorr[t] = pic_corr[t];
          }
          if (pic_corr[t] <= PICcorr[t]) {
            PICcorrSigL[t]++;
          }
          if (pic_corr[t] >= PICcorr[t]) {
            PICcorrSigH[t]++;
          }
          //printf("%f\n", PICcorr[t]);
        }
      }
    }
        
    //calculate variance of standardized contrasts (phylogenetic signal)
    for (k = 0; k<T.ntraits; k++) {
      int count = 0;
      for (node = 0; node < P.nnodes; node++) {
        if (P.noat[node] != 0) {
          c[count] = fabs(contrasts[k][node]);
          count++;
        }
      }
      if (run == 0 && Debug) printf("Node counts:%d\t%d\n",count,P.nnodes-P.termtaxa-1);            
      varC = 0.0;
      for (i = 0;i<count;i++) {
        varC +=  pow(c[i],2);
      }
      varC = varC / (float)(count - 1);         
      if (run == 0) {
        VarC[k] = varC;
      }
      if (varC <= VarC[k]) {
        varCRankLow[k]++;
      }
      if (varC >= VarC[k]) {
        varCRankHigh[k]++;
      }         
    }
  }
  //runs loop just ended

  free_imatrix(RndArr, 0, T.ntaxa, 0, T.ntraits-1);

  if (Output == 0) {
    if (Debug) {
      printf("\nTrait conservatism by node:\n");
      printf("Tr  Nd Ntax  TipMn rank  TipSD rank NdMean rank   NdSD rank Taxon\n");
      for (k = 0; k <= T.ntraits-1; k++) {
        t = ordTrt[k];
        for (node = 0; node < P.nnodes; node++) {
          if (P.noat[node] != 0) {
            printf("%2d %3d ", t+1, node);
            printf("%4d %6.2f ", TipsAtNode[t][node],
                   MeanTipChar[t][node]);
            printf("%4d %6.3f %4d ", TipMeanSigL[t][node],
                   SdTipChar[t][node], TipSdSigL[t][node]);
            printf("%6.2f %4d ", MeanNodeChar[t][node],
                   MeanSigL[t][node]);
            printf("%6.3f %4d ", SdNodeChar[t][node],
                   SdSigL[t][node]);
            printf("%s\n", P.taxon[node]);
          }
        }
      }
    }
    // Print conservatism and PIC results  
    //      printf("\nTrait conservatism - tree wide results:\n");
    //      printf("Tr Ntaxa  MnCont  rank\n");
    //      for (k = 0; k <= T.ntraits-1; k++) {
    //          t = ordTrt[k];
    //          printf("%2d %5d %7.3f %5d\n", t+1, TipsAtNode[t][0], MeanSDTree[t],
    //                  SDTreeSigL[t]);
    //      }
    printf("\nPhylogenetic signal\n");
    printf("Trait\tNTaxa\tVarContr\tVarCn.rankLow\tVarCn.rankHi\n");
    for (k = 0; k <= T.ntraits-1; k++) {
      t = ordTrt[k];
      printf("%s\t%d\t%.3f\t%d\t%d\n", T.trname[t], TipsAtNode[t][0],
             VarC[t],varCRankLow[t],varCRankHigh[t]);
    }
    if (inclContrasts) {
      printf("\nIndependent contrast correlations, traits 2 and higher vs. trait 1:\n");
      if (T.type[ordTrt[0]] == 0) {
        printf("Tr  MnCnt SDCnt  nPos nCont  MnCST");
        printf(" SDCST nPsST nCntST\n");
      } else {
        printf("Tr   PicR  nPos nCont\n");
      }
      for (k = 1; k <= T.ntraits-1; k++) {
        t = ordTrt[k];
        if (T.type[t] != 0) {
          if (T.type[ordTrt[0]] == 0) {
            printf(
                   "%2d %6.3f %5.3f %5.1f %5d %6.3f %5.3f %5.1f  %5d\n",
                   t+1, SumCon[t], SDCon[t], SignCont[t],
                   NSignCont, SumConCons[t], SDConCons[t],
                   SignContCons[t], NSignContCons);
          } else {
            //printf("%2d %6.3f %4d %5.1f %5d %4d\n", t+1,PICcorr[t], PICcorrSigL[t], SignCont[t], nContrast, SignContSigL[t]);  
            printf("%2d %6.3f %5.1f %5d\n", t+1, PICcorr[t],
                   SignCont[t], nContrast);
          }
        }
      }
    }

  }

  if (Output == 1) {
    if ((aotOut == 0) || (aotOut == 1)) {
      //TRAITS OUTPUT
      if (aotOut == 0)
        printf("\nTrait conservatism by node:\n");
      printf("trait\ttrait.name");
      printf("\tnode\tname\tage\tNtaxa\tN.nodes\t");
      printf("Tip.mn\tTmn.rankLow\tTmn.rankHi\t");
      printf("Tip.sd\tTsd.rankLow\tTsd.rankHi\t");
      printf("Node.mn\tNmn.rankLow\tNmn.rankHi\t");
      printf("Nod.sd\tNsd.rankLow\tNsd.rankHi\t");
      if (Debug)
        printf("BLwts\tBL\ttBL\t");
      printf("SSTipsRoot\tSSTips\tpercVarAmongNodes\tpercVarAtNode\tContributionIndex\t");
      printf("SSTipVNodeRoot\tSSTipVNode\tSSAmongNodes\tSSWithinNodes");
      if (Debug)
        printf("\tdepth");
      printf("\n");
      for (k = 0; k <= T.ntraits-1; k++) {
        t = ordTrt[k];
        for (node = 0; node < P.nnodes; node++) {
          if (P.noat[node] != 0) {
            printf("%d\t", t+1);
            printf("%s\t", T.trname[t]);
            printf("%d\t%s", node, P.taxon[node]);
            printf("\t%f", P.age[node]);
            printf("\t%d\t%d\t%.3f\t", TipsAtNode[t][node],
                   P.noat[node], MeanTipChar[t][node]);
            printf("%d\t%d\t%.4f\t%d\t%d\t", TipMeanSigL[t][node],
                   TipMeanSigH[t][node], SdTipChar[t][node],
                   TipSdSigL[t][node], TipSdSigH[t][node]);
            printf("%.3f\t%d\t%d\t", MeanNodeChar[t][node],
                   MeanSigL[t][node], MeanSigH[t][node]);
            printf("%.3f\t%d\t%d\t", SdNodeChar[t][node],
                   SdSigL[t][node], SdSigH[t][node]);
            if (Debug) {
              printf("%f\t%f\t%f\t", BLwts[node], P.bl[node],
                     Tbl[node]);
            }
            printf("%f\t%f\t", SSTipChar[t][0], SSTipChar[t][node]);
            tmpf1 = SSAmongNodes[t][node] / (SSAmongNodes[t][node]
                                             +SSWithinNodes[t][node]);
            tmpf2 = SSTipsVNode[t][node]/SSTipsVNode[t][0];
            printf("%f\t", tmpf1);
            printf("%f\t", tmpf2);
            printf("%f\t", tmpf1*tmpf2);
            printf("%f\t%f\t%f\t%f", SSTipsVNode[t][0],
                   SSTipsVNode[t][node], SSAmongNodes[t][node],
                   SSWithinNodes[t][node]);
            if (Debug)
              printf("\t%d", P.depth[node]);
            printf("\n");
          }
        }
      }
    }
    //This loop prints values of independent contrasts to a 
    //separate block of output file
    if ((aotOut == 0) || (aotOut == 2)) {
      if (inclContrasts)
        if (aotOut == 0)
          printf("\nOutput of independent contrast values by node\n");
      {
        printf("node\tname\tage\tN.nodes");
        for (i = 0; i <= T.ntraits - 1; i++) {
          printf("\tContrast%d", i+1);
        }
        if (T.type[ordTrt[0]] == 0) {
          printf("\tNodeType");
        } else {
          printf("\tContrastSD");
        }
        for (i = 0; i <= T.ntraits - 1; i++) {
          printf("\tLowVal%d\tHiVal%d", i+1, i+1);
        }
        printf("\n");

        for (node = 0; node < P.nnodes; node++) {
          if (nodeStatus[ordTrt[0]][node] == 1) {
            printf("%d\t%s", node, P.taxon[node]);
            printf("\t%f", P.age[node]);
            printf("\t%d", P.noat[node]);
            for (i = 0; i <= T.ntraits - 1; i++) {
              printf("\t%f", Contrast[ordTrt[i]][node]);
            }
            if (T.type[ordTrt[0]] == 0) {
              if (nodeStatusCons[ordTrt[0]][node] == 1) {
                printf("\tST");
              } else {
                printf("\tPT");
              }
            } else {
              printf("\t%f", SDContrast[node]);
            }
            for (i = 0; i <= T.ntraits - 1; i++) {
              printf("\t%f\t%f", DCharLow[ordTrt[i]][node],
                     DCharHi[ordTrt[i]][node]);
            }
            printf("\n");
          }

        }
      }
    }

    if (aotOut == 0) {
      // Loop below prints trait data for Terminals
      if (0) // Set to 1 to print
        {
          if (aotOut == 0)
            printf("\nTrait data for terminal nodes\n");
          printf("node\tname\tnode.depth\ttrait\ttrait.val\n");
          for (node = 0; node < P.nnodes; node++) {
            if (P.noat[node] == 0) {
              for (node = 0; node < P.nnodes; node++)
                for (k = 0; k <= T.ntraits-1; k++) {
                  t = ordTrt[k];
                  printf("%d\t%s\t%d\t%d\t%f\n", node,
                         P.taxon[node], P.depth[node], t+1,
                         MeanNodeChar[t][node]);
                }
            }
          }
        }

      // Output of DOT test

      if (showDOT) {
        printf("\ntrait\tWtdAge\tDOT\tsig\n");
        for (k = 0; k <= T.ntraits-1; k++) {
          t = ordTrt[k];
          printf("%s\t%f\t%f\t%d\n", T.trname[t], wtdContAgeObs[t],
                 DOTObs[t], DOTSig[t]);
        }
        printf("\n");
        // printf("%f\t%d\n",SvSObs[1],SvSSig[1]);
        // printf("%f\t%d\n",DOTObs[1],DOTSig[1]);
        // printf("\n");
      }
    }

    if ((aotOut == 0) || (aotOut == 3)) {

      // Print signal/conservatism and PIC results  
      printf("\nPhylogenetic signal\n");
      printf("Trait\tNTaxa\tVarContr\tVarCn.rankLow\tVarCn.rankHi\n");
      for (k = 0; k <= T.ntraits-1; k++) {
        t = ordTrt[k];
        printf("%s\t%d\t%.3f\t%d\t%d\n", T.trname[t], TipsAtNode[t][0],
               VarC[t],varCRankLow[t],varCRankHigh[t]);
      }
            
    }

    if (inclContrasts && (aotOut == 0)) {
      printf("\nIndependent contrast correlations, traits 2 and higher vs. trait 1:\n");
      if (T.type[ordTrt[0]] == 0) {
        printf("XTrait\tYTrait\tNtaxa\tMnConAll\tSDConAll\tnPosAll\tnContAll\tMnConST\tSDConST\tnPosST\tnContST\n");
      } else {
        printf("XTrait\tYTrait\tNtaxa\tPicR\tnPos\tnCont\n");
      }

      for (k = 1; k <= T.ntraits-1; k++) {
        t = ordTrt[k];
        if (T.type[t] != 0) {
          if (T.type[ordTrt[0]] == 0) {
            printf(
                   "%s\t%s\t%d\t%.3f\t%.3f\t%.1f\t%d\t%.3f\t%.3f\t%.1f\t%d\n",
                   T.trname[ordTrt[0]], T.trname[t],
                   TipsAtNode[t][0], SumCon[t], SDCon[t],
                   SignCont[t], NSignCont, SumConCons[t],
                   SDConCons[t], SignContCons[t], NSignContCons);
          } else {
            printf("%s\t%s\t%d\t%.3f\t%.1f\t%d\n",
                   T.trname[ordTrt[0]], T.trname[t],
                   TipsAtNode[t][0], PICcorr[t], SignCont[t],
                   nContrast);
          }
        }
      }
    }
  }

  if (Output == 2) {
    // BEGIN CAM EDIT
    for (i = 0; i < T.ntraits; i++) {
      Out[i] = P; // all the pointers in Out are the same as those in P
      // careful not to change any arrays in Out, or they will
      // also be changed in P!

      Out[i].arenotes = 1;

      // dimension a new array for Out names:
      Out[i].notes = cmatrix(0, Out[i].nnodes-1, 0, MAXNOTELENGTH+10);
    }

    // for each trait
    for (k = 0; k < T.ntraits; k++) {
      t = ordTrt[k];
      // for each node - not the root node
      for (node = 1; node < P.nnodes; node++) {
        // clear if not yet initialized
        if (Out[t].arenotes == 1)
          strcpy(Out[t].notes[node], "");

        if (TipSdSigL[t][node] <= (int) ((float) RUNS * 0.05)) {
          strcpy(temp, "");
          sprintf(temp, "CONSERVED %d", TipSdSigL[t][node]);
          strcat(Out[t].notes[node], temp);
        }

        if (TipSdSigH[t][node] <= (int) ((float) RUNS * 0.05)) {
          strcpy(temp, "");
          sprintf(temp, "DIVERGENT %d", TipSdSigH[t][node]);
          strcat(Out[t].notes[node], temp);
        }
      }
      // Name tree
      strcpy(Out[t].phyname, "Traits_");
      strcat(Out[t].phyname, T.trname[t]);

      Out[t].arenotes = 2; // notify other functions that notes are initialized
    }
    //printf("]\n\n\n");
    WriteNexus(Out, T.ntraits, ReadSample(SampleFile), 1, T, 1);

    // END CAM EDIT
  }

  free_imatrix(addedAtNode, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(addedAtNodeBin, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(binChar, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(binCharAssign, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(nodeChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(meanNodeCharBin, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(meanTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(nodeStatus, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(nodeStatusCons, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(sdNodeChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(sdTipChar, 0, T.ntraits, 0, P.nnodes);
  free_matrix(sumNodeChar, 0, T.ntraits-1, 0, P.nnodes);
  free_vector(sumSDTree, 0, T.ntraits-1);
  free_matrix(sumTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(sumsqNodeChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(sumsqTipCharG, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(SumsqTipCharG, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(ssBG1, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(ssTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(percVarDivTips, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(PercVarDiv, 0, T.ntraits-1, 0, P.nnodes);
  free_vector(sumsqSDTree, 0, T.ntraits-1);
  free_matrix(sumsqTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_vector(charMean, 0, T.ntraits-1);

  free_vector(depthCorr, 0, T.ntraits-1);
  free_vector(meanSDTree, 0, T.ntraits-1);
  free_vector(sd_SDTree, 0, T.ntraits-1);
  free_vector(sumDepth, 0, T.ntraits-1);
  free_vector(sumNodeSd, 0, T.ntraits-1);
  free_vector(sumprDpNvr, 0, T.ntraits-1);
  free_vector(sumsqDepth, 0, T.ntraits-1);
  free_vector(sumsqNodeSd, 0, T.ntraits-1);
  free_ivector(internalNodes, 0, T.ntraits-1);

  free_matrix(MeanNodeChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(MeanTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(SdNodeChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(SdTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(SSTipChar, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(TipsAtNode, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(MeanSigL, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(SdSigL, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(TipMeanSigL, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(TipSdSigL, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(MeanSigH, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(SdSigH, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(TipMeanSigH, 0, T.ntraits-1, 0, P.nnodes);
  free_imatrix(TipSdSigH, 0, T.ntraits-1, 0, P.nnodes);

  free_matrix(Contrast, 0, T.ntraits-1, 0, P.nnodes);
  free_matrix(contrasts, 0, T.ntraits-1, 0, P.nnodes);  
  free_vector(c, 0, P.nnodes - P.termtaxa - 1);
  free_vector(VarC, 0, T.ntraits-1);
  free_ivector(varCRankHigh, 0, P.nnodes);
  free_ivector(varCRankLow, 0, P.nnodes);
    
}

// ---------- RANDOM ARRAY ----------------------------------------------------

void RandArrayT(int **RndArr, int nTaxa, int run, int nTraits) {
  // Produces a random list from 0 to n-1 inclusive for each trait
  // IF TRAITS ARE KEPT TOGETHER, THEN TEST IS FOR EFFECTS OF PHYLOGENY
  int data[nTaxa];
  int i, j, loop, size, pick, t;

  //printf("%d\n", n);
  for (t = 0; t < nTraits; t++) {
    if (run == 0) {
      for (i = 0; i < nTaxa; i++) {
        RndArr[i][t] = i;
      }
    } else {
      for (i = 0; i < nTaxa; i++) {
        data[i] = i;
      }
      for (loop = 0; loop < nTaxa; loop++) {
        size = nTaxa - loop;
        // pick = ((int) (((float) size * random()) / (RAND_MAX+1.0)));
        pick = (int)((((double) random())/((double) RAND_MAX))
                     *(double) size);
        if (pick == size) {
          pick = pick - 1;
        }
        if (pick > size) {
          printf("RndArr out of bounds, run = %d, loop = %d", run,
                 loop);
        }
        RndArr[loop][t] = data[pick];
        for (j = pick; j < size - 1; j++) {
          data[j] = data[j + 1];
        }
        //printf("%d\t%d\n", t, RndArr[loop][t]);
        //printf("%4d %4d %7.3f\n",loop,RndArr[loop][t],InC.tr[RndArr[loop][t]][t]);
      }
      //if (run == 3)  exit(EXIT_FAILURE);
    }
  }
}
