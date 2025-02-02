/*  Header file for evolve.c                               */

/*  
   Sequence Generator - seq-gen, version 1.3.5
   Copyright (c)1996-2025, Andrew Rambaut
   Institute of Evolutionary Biology, University of Edinburgh			
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://tree.bio.ed.ac.uk/software/seqgen/
   email: a.rambaut@ed.ac.uk
*/


#ifndef _EVOLVE_H_
#define _EVOLVE_H_

#include "tree.h"

#define MAX_RATE_CATS 32

extern int numTaxa, numSites, maxPartitions, numPartitions, fileFormat;

extern double gammaShape, proportionInvariable;
extern int numCats, rateHetero, invariableSites;
extern double catRate[MAX_RATE_CATS];

enum {
	NoRates,
	CodonRates,
	GammaRates,
	DiscreteGammaRates
};

enum {
	PHYLIPFormat,
	RelaxedFormat,
	NEXUSFormat,
	FASTAFormat
};

/* prototypes */

void SetModel(int model);
void SetCategories();
void CreateRates();
void CreateSequences(TTree *tree, int inNumSites);
void EvolveSequences(TTree *tree, int inFromSite, int inNumSites, double scale, char *ancestor);

void WriteSequences(FILE *fv, int treeNo, int datasetNo, TTree **treeSet, int *partitionLengths);

void WriteAncestralSequences(FILE *fv, TTree *tree);

void WriteRates(FILE *fv);

#endif /* _EVOLVE_H_ */
