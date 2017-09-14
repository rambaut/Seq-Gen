/*  
   Sequence Generator - seq-gen, version 1.3.4
   Copyright (c)1996-2017, Andrew Rambaut & Nick Grassly
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "global.h"
#include "tree.h"
#include "evolve.h"
#include "model.h"
#include "gamma.h"
#include "twister.h"

int numTaxa, numSites, maxPartitions, numPartitions, fileFormat;

double gammaShape, proportionInvariable;
int numCats, rateHetero, invariableSites;
double catRate[MAX_RATE_CATS];
double freqRate[MAX_RATE_CATS];

static double *matrix[MAX_RATE_CATS];
static double *vector;

static double *gammaRates;
static short *categories;
static short *invariable;
static double *siteRates;

/* prototypes */

char SetState(double *P);
short IsInvariable();
void SetSequence(char *seq, char *source, int inFromSite, int inNumSites);
void SetNucSequence(char *seq, char *source, int inFromSite, int inNumSites);
void SetAASequence(char *seq, char *source, int inFromSite, int inNumSites);
void RandomSequence(char *seq, int inNumSites);
void MutateSequence(char *seq, int inFromSite, int inNumSites, double len);
void SetSiteRates(int inFromSite, int inNumSites, double inScale);
void EvolveNode(TNode *anc, TNode *des, int inFromSite, int inNumSites, double scale);

void WriteFastaFormat(FILE *fv, TTree **treeSet, int *partitionLengths);
void WritePhylipFormat(FILE *fv, TTree **treeSet, int *partitionLengths);
void WriteNexusFormat(FILE *fv, int treeNo, int datasetNo, TTree **treeSet, int *partitionLengths);

void WriteAncestralSequencesNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des);

/* functions */

void CreateRates()
{
	int i;
	
	if (rateHetero==GammaRates)
		gammaRates=CAllocMem(numSites*sizeof(double), "gammaRates", "CreateRates", 0);
	else if (rateHetero==DiscreteGammaRates)
		categories=CAllocMem(numSites*sizeof(short), "categories", "CreateRates", 0);
		
	if (invariableSites)
		invariable=CAllocMem(numSites*sizeof(short), "invariable", "CreateRates", 0);

	siteRates=CAllocMem(numSites*sizeof(double), "siteRates", "CreateRates", 0);

	for (i = 0; i < MAX_RATE_CATS; i++) {
		matrix[i] = CAllocMem(numStates*numStates*sizeof(double), "matrix", "CreateRates", 0);
	}
	vector = CAllocMem(numStates*sizeof(double), "vector", "CreateRates", 0);
}

void CreateSequences(TTree *tree, int inNumSites)
{
	TNode *P;
	
	P=tree->nodeList;
	while (P!=NULL) {
		if (P->sequence != NULL)
			free(P->sequence);
			
		P->sequence=CAllocMem(inNumSites+1, "sequences", "CreateSequences", 0);
		P=P->next;
	}
}


void SetCategories()
{
	int i;
	double sumRates;
	
	if (rateHetero==CodonRates) {
		sumRates=catRate[0]+catRate[1]+catRate[2];
		if (sumRates!=3.0) {
			catRate[0]*=3.0/sumRates;
			catRate[1]*=3.0/sumRates;
			catRate[2]*=3.0/sumRates;
		}
	} else if (rateHetero==GammaRates) {
		for (i=0; i<numSites; i++)
			gammaRates[i]=rndgamma(gammaShape) / gammaShape;
	} else if (rateHetero==DiscreteGammaRates) {
		DiscreteGamma(freqRate, catRate, gammaShape, gammaShape, numCats, 0);
		for (i=0; i<numSites; i++)
			categories[i]=(int)(rndu()*numCats);
	}

	if (invariableSites) {
		for (i=0; i<numSites; i++)
			invariable[i] = IsInvariable();
	} 
}

/*
Changed "j<numStates" to "j<numStates-1" in the following function to prevent array overflow

20 October 2011 -- Michael Ott, CSIRO (michael.ott@csiro.au)
*/

char SetState(double *P)
{
	char j;
	double r;
	
	r=rndu();
	for (j=0; r>(*P) && j<numStates-1; j++) P++;
	return j;
}

short IsInvariable()
{
	double r;

	r=rndu();
	if (r < proportionInvariable)
		return 1;
	else 
		return 0;
}


void RandomSequence(char *seq, int inNumSites)
{
	int i;
	char *P;
	
	P=seq;
	for (i=0; i<inNumSites; i++) {
		*P=SetState(addFreq);
		P++;
	}
}

void SetSequence(char *seq, char *source, int inFromSite, int inNumSites)
{
	int i, j;
	char *P, *Q;
	
	P=seq;
	Q=source + inFromSite;
	for (i=0; i<inNumSites; i++) {
		for (j=0; j < numStates; j++) {
			if (*Q == stateCharacters[j]) break;
		}
		if (j == numStates) {
			fprintf(stderr, "Bad state in ancestoral sequence\n");
			exit(0);
		}
		*P = j;
		P++;
		Q++;
	}
}

void MutateSequence(char *seq, int inFromSite, int inNumSites, double len)
{
	int i, j, cat;
	double *Q;
	short *R;
	short *S;
	char *P;
	
	P=seq;
	
	/* This is done quite long-windedly to keep the speed up */
	switch (rateHetero) {
		case GammaRates:
			Q=gammaRates + inFromSite;
			
			if (invariableSites) {
				S = invariable + inFromSite;
				for (i=0; i<inNumSites; i++) {
					if (!(*S)) {
						SetVector(vector, *P, (*Q)*len);
						*P=SetState(vector);
					}
					P++;
					Q++;
					S++;
				}
			} else {
				for (i=0; i<inNumSites; i++) {
					SetVector(vector, *P, (*Q)*len);
					*P=SetState(vector);
					P++;
					Q++;
				}
			}
		break;
		case DiscreteGammaRates:
			for (i=0; i<numCats; i++)
				SetMatrix(matrix[i], catRate[i]*len);
			
			R=categories + inFromSite;
			if (invariableSites) {
				S = invariable + inFromSite;
				for (i=0; i<inNumSites; i++) {
					if (!(*S))
						*P=SetState(matrix[*R]+((*P) * numStates));
					P++;
					R++;
					S++;
				}
			} else {
				for (i=0; i<inNumSites; i++) {
					*P=SetState(matrix[*R]+((*P) * numStates));
					P++;
					R++;
				}
			}
		break;
		case CodonRates:
			for (i=0; i<numCats; i++)
				SetMatrix(matrix[i], catRate[i]*len);
			
			j = inFromSite;
			for (i=0; i< inNumSites; i++) {
				cat=j%3;
				*P=SetState(matrix[cat]+((*P) * numStates));
				P++;
				j++;
			}
		break;
		case NoRates:
			SetMatrix(matrix[0], len);
			
			if (invariableSites) {
				S = invariable + inFromSite;
				for (i=0; i<inNumSites; i++) {
					if (!(*S))
						*P=SetState(matrix[0]+((*P) * numStates));
					P++;
					S++;
				}
			} else {
				for (i=0; i<inNumSites; i++) {
					*P=SetState(matrix[0]+((*P) * numStates));
					P++;
				}
			}
		break;
	}
}


void SetSiteRates(int inFromSite, int inNumSites, double inScale)
{
	int i, j, cat;
	double *Q;
	short *R;
	short *S;
	double *P;
	
	P=siteRates + inFromSite;
	
	/* This is done quite long-windedly to keep the speed up */
	switch (rateHetero) {
		case GammaRates:
			Q=gammaRates + inFromSite;
			
			if (invariableSites) {
				S = invariable + inFromSite;
				for (i=0; i<inNumSites; i++) {
					if (!(*S))
						*P = (*Q) * inScale;
					else
						*P = 0.0;
					P++;
					Q++;
					S++;
				}
			} else {
				for (i=0; i<inNumSites; i++) {
					*P = (*Q) * inScale;
					P++;
					Q++;
				}
			}
		break;
		case DiscreteGammaRates:
			R=categories + inFromSite;

			if (invariableSites) {
				S = invariable + inFromSite;
				for (i=0; i<inNumSites; i++) {
					if (!(*S))
						*P = catRate[*R] * inScale;
					else
						*P = 0.0;
					P++;
					R++;
					S++;
				}
			} else {
				for (i=0; i<inNumSites; i++) {
					*P = catRate[*R] * inScale;
					P++;
					R++;
				}
			}
		break;
		case CodonRates:
			j = inFromSite;
			for (i=0; i< inNumSites; i++) {
				cat=j%3;
				*P = catRate[cat] * inScale;
				P++;
				j++;
			}
		break;
		case NoRates:
			if (invariableSites) {
				S = invariable + inFromSite;
				for (i=0; i<inNumSites; i++) {
					if (!(*S))
						*P = inScale;
					else
						*P = 0.0;
					P++;
					S++;
				}
			} else {
				for (i=0; i<inNumSites; i++) {
					*P = inScale;
					P++;
				}
			}
		break;
	}
}


void EvolveNode(TNode *anc, TNode *des, int inFromSite, int inNumSites, double scale)
{
	double len;
	
	len=des->length0*scale;
		
	memcpy(des->sequence, anc->sequence, inNumSites);
	MutateSequence(des->sequence, inFromSite, inNumSites, len);
	
	if (des->tipNo==-1) {
		EvolveNode(des, des->branch1, inFromSite, inNumSites, scale);
		EvolveNode(des, des->branch2, inFromSite, inNumSites, scale);
	}
}


void EvolveSequences(TTree *tree, int inFromSite, int inNumSites, double scale, char *ancestor)
{
	if (ancestor==NULL)
		RandomSequence(tree->root->sequence, inNumSites);
	else
		SetSequence(tree->root->sequence, ancestor, inFromSite, inNumSites);
	
	/* Rescale the branch lengths to accommodate the invariable sites */
	if (invariableSites)
		scale /= (1.0 - proportionInvariable);
	
	SetSiteRates(inFromSite, inNumSites, scale);	
	
	EvolveNode(tree->root, tree->root->branch1, inFromSite, inNumSites, scale);
	EvolveNode(tree->root, tree->root->branch2, inFromSite, inNumSites, scale);
	if (!tree->rooted)
		EvolveNode(tree->root, tree->root->branch0, inFromSite, inNumSites, scale);
}

void WriteSequences(FILE *fv, int treeNo, int datasetNo, TTree **treeSet, int *partitionLengths)
{
	switch (fileFormat) {
		case PHYLIPFormat:
		case RelaxedFormat:
			WritePhylipFormat(fv, treeSet, partitionLengths);
		break;
		case NEXUSFormat:
			WriteNexusFormat(fv, treeNo, datasetNo, treeSet, partitionLengths);
		break;
		case FASTAFormat:
			WriteFastaFormat(fv, treeSet, partitionLengths);
		break;
	}
}

void WritePhylipFormat(FILE *fv, TTree **treeSet, int *partitionLengths)
{
	int i, j, k;
	char *P;

	fprintf(fv, " %d %d\n", numTaxa, numSites);

	for (i=0; i<numTaxa; i++) {
		if (fileFormat == RelaxedFormat)
			fprintf(fv, "%s ", treeSet[0]->names[i]);
		else {
			j=0;
			P=treeSet[0]->names[i];
			while (j<10 && *P) {
				fputc(*P, fv);
				j++;
				P++;
			}
			while (j<10) {
				fputc(' ', fv);
				j++;
			}
		}
		for (k = 0; k < numPartitions; k++) {
			P=treeSet[k]->tips[i]->sequence;
			for (j=0; j<partitionLengths[k]; j++) {
				fputc(stateCharacters[(int)*P], fv);
				P++;
			}
		}
		fputc('\n', fv);
	}
}

void WriteFastaFormat(FILE *fv, TTree **treeSet, int *partitionLengths)
{
	size_t i, j, k;

	for (i=0; i<numTaxa; i++) {
		size_t printedSites = 0;

		fprintf(fv, ">%s\n", treeSet[0]->names[i]);

		for (j = 0; j < numPartitions; j++) {
			char *P;
			P=treeSet[j]->tips[i]->sequence;
			for (k=0; k<partitionLengths[j]; k++) {
				fputc(stateCharacters[(int)*P], fv);
				P++;
				if (++printedSites % 72 == 0) {
					fputc('\n', fv);
				}
			}
		}
		if (printedSites % 72 != 0) {
			fputc('\n', fv);
		}
	}
}

void WriteNexusFormat(FILE *fv, int treeNo, int datasetNo, TTree **treeSet, int *partitionLengths)
{
	int i, j, k, len, maxLen;
	char *P;
	
	if (treeNo > 0 && datasetNo > 0)
		fprintf(fv, "Begin DATA;\t[Tree %d, Dataset %d]\n", treeNo, datasetNo);
	else if (treeNo > 0)
		fprintf(fv, "Begin DATA;\t[Tree %d]\n", treeNo);
	else if (datasetNo > 0)
		fprintf(fv, "Begin DATA;\t[Dataset %d]\n", datasetNo);
	else
		fprintf(fv, "Begin DATA;\n");
		
	fprintf(fv, "\tDimensions NTAX=%d NCHAR=%d;\n", numTaxa, numSites);
	if (isNucModel) {
		fprintf(fv, "\tFormat MISSING=? GAP=- DATATYPE=DNA;\n");
	} else {
		fprintf(fv, "\tFormat MISSING=? GAP=- DATATYPE=PROTEIN;\n");
	}

	fprintf(fv, "\tMatrix\n");

	maxLen = strlen(treeSet[0]->names[0]);
	for (i=1; i<numTaxa; i++) {
		len = strlen(treeSet[0]->names[i]);
		if (len > maxLen)
			maxLen = len;
	}

	for (i=0; i<numTaxa; i++) {
		fprintf(fv, "%s ", treeSet[0]->names[i]);
		len = maxLen - strlen(treeSet[0]->names[i]);
		for (j = 0; j < len; j++) {
			fputc(' ', fv);
		}
		
		for (k = 0; k < numPartitions; k++) {
			P=treeSet[k]->tips[i]->sequence;
			for (j=0; j<partitionLengths[k]; j++) {
				fputc(stateCharacters[(int)*P], fv);
				P++;
			}
		}
		fputc('\n', fv);
	}
	fprintf(fv, "\t;\nEND;\n\n");
}

void WriteAncestralSequences(FILE *fv, TTree *tree)
{
	int j, n;
	char *P;

	if (!tree->rooted)
		n = (2 * numTaxa) - 3;
	else
		n = (2 * numTaxa) - 2;
	fprintf(fv, " %d %d\n", n, numSites);

	n = numTaxa + 1;
	
	fprintf(fv, "%d\t", n);
	P=tree->root->sequence;
	for (j=0; j<numSites; j++) {
		fputc(stateCharacters[(int)*P], fv);
		P++;
	}
	fputc('\n', fv);
	
	if (!tree->rooted)
		WriteAncestralSequencesNode(fv, tree, &n, tree->root->branch0);
	WriteAncestralSequencesNode(fv, tree, &n, tree->root->branch1);
	WriteAncestralSequencesNode(fv, tree, &n, tree->root->branch2);
}

void WriteAncestralSequencesNode(FILE *fv, TTree *tree, int *nodeNo, TNode *des)
{
	int j;
	char *P = des->sequence;
	
	if (des->tipNo==-1) {
		(*nodeNo)++;
		
		fprintf(fv, "%d\t", *nodeNo);
		for (j=0; j<numSites; j++) {
			fputc(stateCharacters[(int)*P], fv);
			P++;
		}
		fputc('\n', fv);
		
		WriteAncestralSequencesNode(fv, tree, nodeNo, des->branch1);
		WriteAncestralSequencesNode(fv, tree, nodeNo, des->branch2);
	} else {
		fprintf(fv, "%s\t", tree->names[des->tipNo]);
		for (j=0; j<numSites; j++) {
			fputc(stateCharacters[(int)*P], fv);
			P++;
		}
		fputc('\n', fv);
	}

}


void WriteRates(FILE *fv)
{
	int i;
	
	fprintf(fv, "Relative rates for each site:\n");
	
	for (i=0; i<numSites; i++) {
		fprintf(fv, "%d\t%lf\n", i+1, siteRates[i]);
	}
	fprintf(fv, "\n");
}


