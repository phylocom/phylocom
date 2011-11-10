/*  ecovolve.h  */

// This program is free software; you can redistribute it and/or
// modify it under the terms of the BSD 2-Clause License
//
//   http://www.opensource.org/licenses/bsd-license.php
//
// Copyright (c) 2009-2011, Campbell Webb <cwebb@oeb.harvard.edu>
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

// Routine:
// new time
//   {
//      determine competetive environ
//      (begin real time)
//      evolve chars on extant lineages
//      (end real time)
//      speciate; new spp get old state; make mid-time marker in t+1 on.
//      extinct; active off (newly speciated clades cannot go extinct)
//        else: pass mid-time lineage markers to future 
//   }

// I.e. for time x lineage is before node


/* INCLUDE HEADERS ------------------------------------------------------- */

#include <stdio.h>
#include <math.h>  // for sqrt()
#include <stdlib.h>
#include "phylocom.h"

/* DEFINITIONS FOR MAIN PROGRAM ------------------------------------------ */

#define OUTMODE 3  // 0 = debug, 1 = .phy out, 2 = lineage-thr-time
#define P_SPECIATE 0.05
#define P_EXTINCT 0.01
//#define COMPETE 0
//#define MAXNODES 50000 // now defined in phylocom
#define MAXTIME 100


// for fy2new
//#define MAXTAXONLENGTH 100
//#define MAXBLLENGTH 15
//#define MAXNOTELENGTH 100

/* FUNCTION DECLARATION -------------------------------------------------- */

void Speciate();
void Extinct();
void CharChange();
void Output();
void Compete();
void MakePhylo();
void WriteTraits();
void DummySample();
struct shift MakeChange();
//struct phylo Prune();  // maybe move to phylocom.h?
float Balance();

/* GLOBAL VARIABLES ------------------------------------------------------ */

int OUT_MODE;
int MAX_TIME;
int PRUNED;
int NCHAR;
int TAPER = 0;
int TAPERFACT = 1;
int COMPETE = 0;

int Lineage; // incremental counter of higher Lineage no
int Node;    // ditto for node
int Name;    // ditto for name

// tree tracing:
int Active_l[MAXNODES]; // active marker during tree gen
int Living_lt[MAXNODES][MAXTIME+1]; // tracer for tree history
int Censor_lt[MAXNODES][MAXTIME+1];     // tracer for cenored tree
int EcoDist_l[MAXNODES];

// lineage states
int LineageUp_l[MAXNODES];          
int NodeUp_l[MAXNODES];
float ***Char_ltc;
int Extant_l[MAXNODES];

// node states:
int Time_n[MAXNODES];
int NodeUp_n[MAXNODES];
int BlUp_n[MAXNODES];
int Name_n[MAXNODES];
float **Char_nc;
int Terminal_n[MAXNODES];

float *CharMin, *CharMax;
float *CharMean, *CharSumXSqr;
float *CharCMin, *CharCMax;
float *CharCMean, *CharCSumXSqr;
int LttZ;
int LttC;

phylo OutTree;
// phylo PrunedTree;

typedef struct shift {
  int n;
  int *x;
} shift;

shift charch;

