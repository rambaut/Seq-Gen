/*  
   Sequence Generator - seq-gen, version 1.3.2
   Copyright (c)1996-2004, Andrew Rambaut & Nick Grassly
   Department of Zoology, University of Oxford			
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
   http://evolve.zoo.ox.ac.uk/software/Seq-Gen/
   email: andrew.rambaut@zoo.ox.ac.uk
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "evolve.h"
#include "model.h"
#include "nucmodels.h"
#include "aamodels.h"

char *modelNames[numModels]={
	"HKY",
	"F84",
	"GTR",
	"JTT",
	"WAG",
	"PAM",
	"BLOSUM",
	"MTREV",
	"CPREV",
	"GENERAL"
};

char *modelTitles[numModels]={
	"HKY: Hasegawa, Kishino & Yano (1985)",
	"F84: Felsenstein (1984)",
	"GTR: General time reversible (nucleotides)",
	"JTT: Jones, Taylor & Thornton (1992) CABIOS  8:275-282\n             DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
	"WAG: Whelan & Goldman (2001) Mol Biol Evol 18:691Ð699",
	"PAM: Dayhoff, Schwartz & Orcutt (1978)\n             DCMut version Kosiol & Goldman (2004) <http://www.ebi.ac.uk/goldman-srv/dayhoff/>",
	"BLOSUM62: Henikoff & Henikoff (1992) PNAS USA 89:10915-10919",
	"MTREV24: Adachi & Hasegawa (1996) J Mol Evol 42:459-468",
	"CPREV45: Adachi et al. (2000) J Mol Evol 50:348-358",
	"GENERAL: General time reversible (amino acids)"
};

int model, numStates, isNucModel, userFreqs, equalFreqs;

char *stateCharacters;

double *freq, *addFreq;

void SetupFrequencies();
void SetupMatrices();

void InitialiseEigen();
void elmhes(double** a, int* ordr, int n);
void mcdiv(double ar, double ai, double br, double bi);
void hqr2(int n, int low, int hgh, double** h, double** zz, double* wr, double* wi);
void eltran(double** a, double** zz, int* ordr, int n);
void luinverse(double** inmat, double** imtrx, int size);

/*************************************/
void SetModel(int theModel)
{	
	int i;
	
	model=theModel;

	if (isNucModel) {
		numStates = NUM_NUC;
		
		SetNucModel(theModel);

		freq = nucFreq;
		addFreq = nucAddFreq;
		stateCharacters = nucleotides;
	} else {
		numStates = NUM_AA;
		
		SetAAModel(theModel - numNucModels);

		freq = aaFreq;
		addFreq = aaAddFreq;
		stateCharacters = aminoAcids;
	}
	
	addFreq[0]=freq[0];
	for (i = 1; i < numStates; i++) {
		addFreq[i] = addFreq[i-1] + freq[i];
	}
}


void SetMatrix(double *matrix, double len)
{
	if (isNucModel) {
		SetNucMatrix(matrix, len);
	} else {
		SetAAMatrix(matrix, len);
	}
}

void SetVector(double *vector, short state, double len)
{
	if (isNucModel) {
		SetNucVector(vector, state, len);
	} else {
		SetAAVector(vector, state, len);
	}
}

