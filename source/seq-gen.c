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
#include <ctype.h>
#include <time.h>
#include <math.h>

#include "global.h"
#include "treefile.h"
#include "evolve.h"
#include "model.h"
#include "nucmodels.h"
#include "aamodels.h"
#include "progress.h"
#include "twister.h"

#define PROGRAM_NAME "seq-gen"
#define VERSION_NUMBER "Version 1.3.4"

int treeFile, textFile, numDatasets, numTrees;
int scaleTrees, scaleBranches, ancestorSeq, writeAncestors, writeRates;
int *partitionLengths;
double *partitionRates;
double treeScale, branchScale;
char treeFileName[256];
char textFileName[256];

int hasAlignment, numSequences, numAlignmentSites;
char **names;
char **sequences;

FILE *tree_fv;

/* prototypes */

static void PrintTitle();
static void PrintUsage();
static void PrintVerbose(FILE *fv);
static void ReadParams();
static void ReadFileParams();
static void AllocateMemory();
static void ReadFile();
static int OpenTreeFile();
 
/* functions */

static void PrintTitle()
{
	fprintf(stderr, "Sequence Generator - %s\n", PROGRAM_NAME);
	fprintf(stderr, "%s\n", VERSION_NUMBER);
	fprintf(stderr, "(c) Copyright, 1996-2017 Andrew Rambaut and Nick Grassly\n");
	fprintf(stderr, "Institute of Evolutionary Biology, University of Edinburgh\n\n");
	fprintf(stderr, "Originally developed at:\n");
	fprintf(stderr, "Department of Zoology, University of Oxford\n\n");
}

static void PrintUsage()
{ 
	fprintf(stderr, "Usage: seq-gen [-m MODEL] [-l #] [-n #] [-p #] [-s # | -d #] [-k #]\n");
	fprintf(stderr, "               [-c #1 #2 #3 | -a # [-g #]] [-i #] [-f e | #] [-t # | -r #]\n");
	fprintf(stderr, "               [-z #] [-o[p][r][n]] [-w[a][r]] [-x NAME] [-q] [-h] [treefile]\n");
	fprintf(stderr, "  -l: # = sequence length [default = 1000].\n");
	fprintf(stderr, "  -n: # = simulated datasets per tree [default = 1].\n");
	fprintf(stderr, "  -p: # = number of partitions (and trees) per sequence [default = 1].\n");
	fprintf(stderr, "  -s: # = branch length scaling factor [default = 1.0].\n");
	fprintf(stderr, "  -d: # = total tree scale [default = use branch lengths].\n");
	fprintf(stderr, "  -k: # = use sequence k as ancestral (needs alignment) [default = random].\n");

	fprintf(stderr, "\n Substitution model options:\n");
	fprintf(stderr, "  -m: MODEL = HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, MTREV, CPREV45, MTART, LG, GENERAL\n");
	fprintf(stderr, "      HKY, F84 & GTR are for nucleotides the rest are for amino acids\n");
	fprintf(stderr, "  -a: # = shape (alpha) for gamma rate heterogeneity [default = none].\n");
	fprintf(stderr, "  -g: # = number of gamma rate categories [default = continuous].\n");
	fprintf(stderr, "  -i: # = proportion of invariable sites [default = 0.0].\n");

	fprintf(stderr, "\n Nucleotide model specific options:\n");
	fprintf(stderr, "  -c: #1 #2 #3 = rates for codon position heterogeneity [default = none].\n");
	fprintf(stderr, "  -t: # = transition-transversion ratio [default = equal rate].\n");
	fprintf(stderr, "  -r: #1 #2 #3 #4 #5 #6= general rate matrix [default = all 1.0].\n");
	fprintf(stderr, "  -f: #A #C #G #T = nucleotide frequencies [default = all equal].\n");
	
	fprintf(stderr, "\n Amino Acid model specific options:\n");
	fprintf(stderr, "      specify using the order ARNDCQEGHILKMFPSTWYV\n");
	fprintf(stderr, "  -r: #1 .. #190 = general rate matrix [default = all 1.0].\n");
	fprintf(stderr, "  -f: #1 .. #20 = amino acid frequencies e=equal [default = matrix freqs].\n");
	
	fprintf(stderr, "\n Miscellaneous options:\n");
	fprintf(stderr, "  -z: # = seed for random number generator [default = system generated].\n");
	fprintf(stderr, "  -o: Output file format [default = PHYLIP]\n");
	fprintf(stderr, "    p PHYLIP format\n");
	fprintf(stderr, "    r relaxed PHYLIP format\n");
	fprintf(stderr, "    n NEXUS format\n");
	fprintf(stderr, "    f FASTA format\n");
	fprintf(stderr, "  -w: Write additional information [default = none]\n");
	fprintf(stderr, "    a Write ancestral sequences for each node\n");
	fprintf(stderr, "    r Write rate for each site\n");
	fprintf(stderr, "  -x: NAME = a text file to insert after every dataset [default = none].\n");
	fprintf(stderr, "  -h: Give this help message\n");
	fprintf(stderr, "  -q: Quiet\n");
	fprintf(stderr, "  treefile: name of tree file [default = trees on stdin]\n\n");
}


void ReadParams(int argc, char **argv)
{
	int i, j, k;
	char ch, *P, st[4];
	int modelTwoArgs = 0;
	
	model=NONE;

	scaleTrees=0;
	treeScale=0.0;
	scaleBranches=0;
	branchScale=0.0;
	
	maxPartitions=1;
	numPartitions=1;
	
	userSeed = 0;
	
	numCats=1;
	rateHetero=NoRates;
	catRate[0]=1.0;
	gammaShape=1.0;
	
	invariableSites=0;
	proportionInvariable = 0.0;
		
	equalFreqs = 1;
	equalTstv = 1;
	tstv=0.50002;
	
	for (i = 0; i < NUM_AA_REL_RATES; i++) {
		aaRelativeRate[i] = 1.0;
	}
	
	for (i = 0; i < NUM_AA; i++) {
		aaFreq[i] = 1.0;
	}
	
	aaFreqSet = 0;
		
	numSites=-1;
	numDatasets=1;
	
	ancestorSeq=0;
	
	writeAncestors=0;
	writeRates=0;
	
 	verbose=1;
	fileFormat = PHYLIPFormat;
	quiet=0;

	treeFile=0;
	textFile=0;
		
	for (i=1; i<argc; i++) {
		P=argv[i];
		if (*P=='-') {
			P++;
			ch=toupper(*P);
			P++;
			switch (ch) {
				case 'H':
					PrintTitle();
					PrintUsage();
					exit(0);
				break;
				case 'M':
					k = i;
					if (GetStrParam(argc, argv, &i, P, st, 3)) {
						fprintf(stderr, "Bad (or missing) Model Code: %s\n\n", argv[i]);
						exit(1);
					}
					P=st;
					if (i > k) {
						modelTwoArgs = 1;
					}
								
					model=-1;
					for (j=F84; j<numModels; j++) {
						if (strncmp(P, modelNames[j], 3)==0) {
							model=j;
							if (model <= GTR) {
								isNucModel = 1;
								numStates = 4;
							} else {
								isNucModel = 0;
								numStates = 20;
							}
						} else if (strncmp(P, "REV", 3)==0) {
							model=GTR;
							isNucModel = 1;
							numStates = 4;
						}

					}
					if (model==-1) {
						fprintf(stderr, "Unknown Model: %s\n\n", argv[i]);
						exit(1);
					}

				break;
			}
		}
	}

	if (model==NONE) {
		fprintf(stderr, "No model has been specified (use the -m option)\n\n");
		PrintUsage();
		exit(1);
	}

	for (i=1; i<argc; i++) {
		P=argv[i];
		if (*P!='-') {
			if (treeFile) {
				fprintf(stderr, "Illegal command parameter: %s\n\n", argv[i]);
				PrintUsage();
				exit(1);
			}
			treeFile=1;
			strcpy(treeFileName, argv[i]);
		} else if (*P=='-' && toupper(*(P+1))=='X') {
			P++; P++;
			if (*P=='\0') {
				i++;
				P = argv[i];
			}
			textFile=1;
			strcpy(textFileName, P);
		} else {
			P++;
			ch=toupper(*P);
			P++;
			switch (ch) {
				case 'H':
					// already delt with
				break;
				case 'M':
					// already delt with
					if (modelTwoArgs) {
						// the model took two arguments so skip the second one.
						i++;
					}
				break;
				case 'L':
					if (GetIntParams(argc, argv, &i, P, 1, &numSites) || numSites<1) {
						fprintf(stderr, "Bad (or missing) sequence length: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'N':
					if (GetIntParams(argc, argv, &i, P, 1, &numDatasets) || numDatasets<1) {
						fprintf(stderr, "Bad (or missing) number of datasets: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'P':
					if (GetIntParams(argc, argv, &i, P, 1, &maxPartitions) || maxPartitions < 1) {
						fprintf(stderr, "Bad number of partitions: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'C':
					if (!isNucModel) {
						fprintf(stderr, "You can only have codon rates when using nucleotide models\n\n");
						exit(1);
					}
					if (rateHetero==GammaRates) {
						fprintf(stderr, "You can only have codon rates or gamma rates not both\n\n");
						exit(1);
					}
					numCats=3;
					rateHetero=CodonRates;
					if (GetDoubleParams(argc, argv, &i, P, 3, catRate) ||
						catRate[0] <= 0.0 || catRate[1] <= 0.0 || catRate[2] <= 0.0 ) {
						fprintf(stderr, "Bad Category Rates: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'I':
					if (GetDoubleParams(argc, argv, &i, P, 1, &proportionInvariable) || 
							proportionInvariable < 0.0 || proportionInvariable >= 1.0) {
						fprintf(stderr, "Bad Proportion of Invariable Sites: %s\n\n", argv[i]);
						exit(1);
					}
					invariableSites = 1;
				break;
				case 'A':
					if (rateHetero==CodonRates) {
						fprintf(stderr, "You can only have codon rates or gamma rates not both\n\n");
						exit(1);
					}
					
					if (rateHetero==NoRates)
						rateHetero=GammaRates;
					if (GetDoubleParams(argc, argv, &i, P, 1, &gammaShape) || gammaShape<=0.0) {
						fprintf(stderr, "Bad Gamma Shape: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'G':
					if (rateHetero==CodonRates) {
						fprintf(stderr, "You can only have codon rates or gamma rates not both\n\n");
						exit(1);
					}
					
					rateHetero=DiscreteGammaRates;
					if (GetIntParams(argc, argv, &i, P, 1, &numCats) || numCats<2 || numCats>MAX_RATE_CATS) {
						fprintf(stderr, "Bad number of Gamma Categories: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'F':
					if (isNucModel) {
						if (toupper(*P)=='E'){
							/* do nothing - equal freqs is default for nucleotides */
						} else {
							equalFreqs = 0;
							if (GetDoubleParams(argc, argv, &i, P, NUM_NUC, nucFreq)) {
								fprintf(stderr, "Bad Nucleotide Frequencies: %s\n\n", argv[i]);
								exit(1);
							}
						}
					} else {
						aaFreqSet = 1;
						if (toupper(*P)=='E'){
							equalFreqs = 1;
							for(j=0;j<NUM_AA;j++) {
								aaFreq[j]=0.05;
							}
						} else {
							equalFreqs = 0;
							if (GetDoubleParams(argc, argv, &i, P, NUM_AA, aaFreq)) {
								fprintf(stderr, "Bad Amino Acid Frequencies: %s\n\n", argv[i]);
								exit(1);
							}
						}
					}
				break;
				case 'T':
					if (model != HKY && model != F84) {
						fprintf(stderr, "You can only have a transition/transversion ratio when using HKY or F84 models\n\n");
						exit(1);
					}
					equalTstv = 0;
					if (GetDoubleParams(argc, argv, &i, P, 1, &tstv)) {
						fprintf(stderr, "Bad Transition-Transversion Ratio: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'R':
					if (model == GTR) {
						if (GetDoubleParams(argc, argv, &i, P, NUM_NUC_REL_RATES, nucRelativeRates)) {
							fprintf(stderr, "Bad General Nucleotide Rate Matrix: %s\n\n", argv[i]);
							exit(1);
						}
						if (nucRelativeRates[NUM_NUC_REL_RATES - 1]!=1.0) {
							for (j=0; j < NUM_NUC_REL_RATES - 1; j++) 
								nucRelativeRates[j] /= nucRelativeRates[NUM_NUC_REL_RATES - 1];
							nucRelativeRates[NUM_NUC_REL_RATES - 1] = 1.0;
						}
					} else if ( model == GENERAL) {
						if (GetDoubleParams(argc, argv, &i, P, NUM_AA_REL_RATES, aaRelativeRate)) {
							fprintf(stderr, "Bad General Amino Acid Rate Matrix: %s\n\n", argv[i]);
							exit(1);
						}
					} else {
						fprintf(stderr, "You can only have a general rate matrix when using GTR or GENERAL models\n\n");
						exit(1);
					}
				break;
				case 'D':
					scaleTrees=1;
					if (GetDoubleParams(argc, argv, &i, P, 1, &treeScale) || treeScale<=0.0) {
						fprintf(stderr, "Bad Total Tree Scale: %s\n\n", argv[i]);
						exit(1);
					}
					if (scaleBranches) {
						fprintf(stderr, "You can't specify both the -d and -s options\n\n");
						exit(1);
					}
				break;
				case 'S':
					scaleBranches=1;
					if (GetDoubleParams(argc, argv, &i, P, 1, &branchScale) || branchScale<=0.0) {
						fprintf(stderr, "Bad Branch Length Scale: %s\n\n", argv[i]);
						exit(1);
					}
					if (scaleTrees) {
						fprintf(stderr, "You can't specify both the -d and -s options\n\n");
						exit(1);
					}
				break;
				case 'K':
					if (GetIntParams(argc, argv, &i, P, 1, &ancestorSeq) || ancestorSeq<1) {
						fprintf(stderr, "Bad ancestral sequence number: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'Z':
					userSeed = 1;
					if (GetUnsignedLongParams(argc, argv, &i, P, 1, &randomSeed)) {
						fprintf(stderr, "Bad random number generator seed: %s\n\n", argv[i]);
						exit(1);
					}
				break;
				case 'O':
					switch (toupper(*P)) {
						case 'P': fileFormat=PHYLIPFormat; break;
						case 'R': fileFormat=RelaxedFormat; break;
						case 'N': fileFormat=NEXUSFormat; break;
						case 'F': fileFormat=FASTAFormat; break;
						default:					
							fprintf(stderr, "Unknown output format: %s\n\n", argv[i]);
							PrintUsage();
							exit(1);
					}
				break;
				case 'W':
					switch (toupper(*P)) {
						case 'A': writeAncestors=1; break;
						case 'R': writeRates=1; break;
						default:					
							fprintf(stderr, "Unknown write mode: %s\n\n", argv[i]);
							PrintUsage();
							exit(1);
					}
				break;
				case 'Q':
					quiet=1;
				break;
				default:
					fprintf(stderr, "Illegal command parameter: %s\n\n", argv[i]);
					PrintUsage();
					exit(1);
				break;
			}
		}
	}
}

void PrintVerbose(FILE *fv)
{
	int i;
	
	if (numStates == 4)
		fprintf(fv, "Simulations of %d taxa, %d nucleotides\n", numTaxa, numSites);
	else 
		fprintf(fv, "Simulations of %d taxa, %d amino acids\n", numTaxa, numSites);
	fprintf(fv, "  for %d tree(s) with %d dataset(s) per tree\n", numTrees, numDatasets);
	if (numPartitions > 1) {
		fprintf(fv, "  and %d partitions (and trees) per dataset\n", numPartitions);
		fprintf(fv, "    Partition  No. Sites  Relative Rate\n");
		for (i = 0; i < numPartitions; i++) 
			fprintf(fv, "    %4d       %7d    %lf\n", i+1, partitionLengths[i], partitionRates[i]);
	}

	fputc('\n', fv);
	
	if (hasAlignment) {
		fprintf(fv, "Alignment read: numSequences = %d, numAlignmentSites = %d\n", numSequences, numAlignmentSites);
		if (ancestorSeq > 0) {
			fprintf(fv, "Using sequence %d as the ancestral sequence\n", ancestorSeq);
		}
		fputc('\n', fv);
	}
	
	if (scaleTrees) {
		fprintf(fv, "Branch lengths of trees scaled so that tree is %G from root to tip\n\n", treeScale);
	} else if (scaleBranches) {
		fprintf(fv, "Branch lengths of trees multiplied by %G\n\n", branchScale);
	} else {
		fprintf(fv, "Branch lengths assumed to be number of substitutions per site\n\n");
	}
	if (rateHetero==CodonRates) {
		fprintf(fv, "Codon position rate heterogeneity:\n");
		fprintf(fv, "    rates = 1:%f 2:%f 3:%f\n", catRate[0], catRate[1], catRate[2]);
	} else if (rateHetero==GammaRates) {
		fprintf(fv, "Continuous gamma rate heterogeneity:\n");
		fprintf(fv, "    shape = %f\n", gammaShape);
	} else if (rateHetero==DiscreteGammaRates) {
		fprintf(fv, "Discrete gamma rate heterogeneity:\n");
		fprintf(fv, "    shape = %f, %d categories\n", gammaShape, numCats);
	} else
		fprintf(fv, "Rate homogeneity of sites.\n");
	if (invariableSites) {
		fprintf(fv, "Invariable sites model:\n");
		fprintf(fv, "    proportion invariable = %f\n", proportionInvariable);
	}
	fprintf(fv, "Model = %s\n", modelTitles[model]);
	if (isNucModel) {
		if (equalTstv) {
			fprintf(fv, "  Rate of transitions and transversions equal:\n");
		}
		if (model==F84) {
			fprintf(fv, "  transition/transversion ratio = %G (K=%G)\n", tstv, kappa);
		} else if (model==HKY) {
			fprintf(fv, "  transition/transversion ratio = %G (kappa=%G)\n", tstv, kappa);
		} else if (model==GTR) {
			fprintf(fv, "  rate matrix = gamma1:%7.4f alpha1:%7.4f  beta1:%7.4f\n", nucRelativeRates[0], nucRelativeRates[1], nucRelativeRates[2]);
			fprintf(fv, "                                beta2:%7.4f alpha2:%7.4f\n", nucRelativeRates[3], nucRelativeRates[4]);
			fprintf(fv, "                                              gamma2: %7.4f\n", nucRelativeRates[5]);
		}

		if (equalFreqs) {
			fprintf(fv, "  with nucleotide frequencies equal.\n");
		} else {
			fprintf(fv, "  with nucleotide frequencies specified as:\n");
			fprintf(fv, "  A=%G C=%G G=%G T=%G\n\n", freq[A], freq[C], freq[G], freq[T]);
		}
	} else {
		if (aaFreqSet) {
			if (equalFreqs) {
				fprintf(fv, "  with amino acid frequencies equal.\n\n");
			} else {
				fprintf(fv, "  with amino acid frequencies specified as:\n");
				fprintf(fv, "  ");
				for (i = 0; i < NUM_AA; i++) {
					fprintf(fv, " %c=%G", aminoAcids[i], freq[i]);
				}
				fprintf(fv, "\n\n");
			}
		}
	}
}

void ReadFileParams()
{
	char ch, st[256];
	
	hasAlignment=0;
	
	ch=fgetc(tree_fv);
	while (!feof(tree_fv) && isspace(ch)) {
		ch=fgetc(tree_fv);
	}
		
	ungetc(ch, tree_fv);

	if (ch!='(' && isdigit(ch)) {
		fgets(st, 255, tree_fv);
		if ( sscanf( st, " %d %d", &numSequences, &numAlignmentSites)!=2 ) {
			fprintf(stderr, "Unable to read parameters from standard input\n");
			exit(2);
		}

		hasAlignment=1;
		
//		fprintf(stderr, "%d sequences, %d sites\n", numSequences, numAlignmentSites);
	}		
}

void AllocateMemory()
{
	int i;
	
	names=(char **)AllocMem(sizeof(char *)*numSequences, "names", "AllocateMemory", 0);
	sequences=(char **)AllocMem(sizeof(char *)*numSequences, "sequences", "AllocateMemory", 0);
	for (i=0; i<numSequences; i++) {
		names[i]=(char *)AllocMem(sizeof(char)*(MAX_NAME_LEN+1), "names[]", "AllocateMemory", 0);
		sequences[i]=(char *)AllocMem(sizeof(char)*numAlignmentSites, "sequences[]", "AllocateMemory", 0);
	}
}


void ReadFile()
{
	int n, b, i;
	char ch;
		
	n=0;
	do {
		ch=fgetc(tree_fv);
		while ( !feof(tree_fv) && isspace(ch)) {
			ch=fgetc(tree_fv);
		}
			
		if ( feof(tree_fv) ) {
			fprintf(stderr, "Unexpected end of file on standard input\n"); 
			exit(2);
		}
			
		i=0;
		while ( i<MAX_NAME_LEN && !feof(tree_fv) && !isspace(ch) ) {
			names[n][i]=ch;
			ch=fgetc(tree_fv);
			i++;
		}
		names[n][i]='\0';
		if (i==0) {
			fprintf(stderr, "Name missing for species %d\n", n+1);
			exit(2);
		}
		while (!feof(tree_fv) && isspace(ch)) {
			ch=fgetc(tree_fv);
		}
		
		if ( feof(tree_fv) ) {
			fprintf(stderr, "Unexpected end of file on standard input\n");
			exit(2);
		}
		
		b=0;
		while ( !feof(tree_fv) && b<numAlignmentSites) {
			if ( !isspace(ch) ) {
				sequences[n][b]=ch;
				b++;
			}
			ch=toupper(fgetc(tree_fv));
		}

		if ( b<numAlignmentSites ) {
			fprintf(stderr, "Unexpected end of file on standard input\n");
			exit(2);
		}
		
//		fprintf(stderr, "%d: %s, bases read: %d, %s\n", n+1, names[n], b, sequences[n]); 

		n++;
		
		if ( n<numSequences && feof(tree_fv) ) {
			fprintf(stderr, "Too few sequences in input file\n");
			exit(2);
		}
	} while ( n<numSequences );
	
	while (!feof(tree_fv) && isspace(ch)) {
		ch=fgetc(tree_fv);
	}
	ungetc(ch, tree_fv);
}

int OpenTreeFile()
{
	char st[256];
	int n;
		
	if (treeFile) {
		if ( (tree_fv=fopen(treeFileName, "rt"))==NULL ) {
			fprintf(stderr, "Error opening tree file: '%s'\n", treeFileName);
			exit(3);
		}
		n=CountTrees(tree_fv);
	} else {
		tree_fv=stdin;
		if (hasAlignment) {
			fgets(st, 255, stdin);
			if ( sscanf(st, " %d ", &n)!=1 ) {
				fprintf(stderr, "Tree is missing from end of sequence file\n");
				exit(3);
			}
		} else
			n=CountTrees(stdin);
	}
	
	return n;
}

int main(int argc, char **argv)
{
	int i, j, k, treeNo, sumLength;
	char ch;
	TTree **treeSet;
	FILE *text_fv;
	clock_t totalStart;
	double totalSecs, scale, sum;
	char *ancestor;

	totalStart = clock();

	ReadParams(argc, argv);

	if (rateHetero == CodonRates && invariableSites) {
		fprintf(stderr, "Invariable sites model cannot be used with codon rate heterogeneity.\n");
		exit(4);
	}

	if (writeAncestors && fileFormat == NEXUSFormat) {
		fprintf(stderr, "Warning - When writing ancestral sequences, relaxed PHYLIP format is used.\n");
	}

	if (writeAncestors && maxPartitions > 1) {
		fprintf(stderr, "Writing ancestral sequences can only be used for a single partition.\n");
		exit(4);
	}
			
	if (!userSeed)
		randomSeed = CreateSeed();
		
	SetSeed(randomSeed);

	if (!quiet)
 		PrintTitle();
	
	numTrees = OpenTreeFile();

	/* if (!treeFile) { */
		ReadFileParams();
	/*} */


	if ((ancestorSeq>0 && !hasAlignment) || ancestorSeq>numSequences) {
		fprintf(stderr, "Bad ancestral sequence number: %d (%d sequences loaded)\n", ancestorSeq, numSequences);
		exit(4);
	}
	
	if (textFile) {
		if ( (text_fv=fopen(textFileName, "rt"))==NULL ) {
			fprintf(stderr, "Error opening text file for insertion into output: '%s'\n", textFileName);
			exit(4);
		}
	}

	ancestor=NULL;
	if (hasAlignment) {
		AllocateMemory();	
		ReadFile();
		
		if (numSites<0)
			numSites=numAlignmentSites;		
			
		if (ancestorSeq>0) {
			if (numSites!=numAlignmentSites) {
				fprintf(stderr, "Ancestral sequence is of a different length to the simulated sequences (%d)\n", numAlignmentSites);
				exit(4);
			}
			ancestor=sequences[ancestorSeq-1];
		}
	} else if (numSites<0)
		numSites=1000;
	
	SetModel(model);
	
	numTaxa=-1;
	scale=1.0;
	
	treeSet = (TTree **)malloc(sizeof(TTree **) * maxPartitions);
	if (treeSet==NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(5);
	}
	
	partitionLengths = (int *)malloc(sizeof(int) * maxPartitions);
	if (partitionLengths==NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(5);
	}
	
	partitionRates = (double *)malloc(sizeof(double) * maxPartitions);
	if (partitionRates==NULL) {
		fprintf(stderr, "Out of memory\n");
		exit(5);
	}
	
	for (i = 0; i < maxPartitions; i++) {
		if ((treeSet[i]=NewTree())==NULL) {
			fprintf(stderr, "Out of memory\n");
			exit(5);
		}
	}
					
	CreateRates();
	
	treeNo=0;
	do {
		partitionLengths[0] = -1;
		ReadTree(tree_fv, treeSet[0], treeNo+1, 0, NULL, &partitionLengths[0], &partitionRates[0]);

		if (treeNo==0) {
			numTaxa=treeSet[0]->numTips;
			
			if (!quiet)
				fprintf(stderr, "Random number generator seed: %ld\n\n", randomSeed);
				
			if (fileFormat == NEXUSFormat) {
				fprintf(stdout, "#NEXUS\n");
				fprintf(stdout, "[\nGenerated by %s %s\n\n", PROGRAM_NAME, VERSION_NUMBER);
				PrintVerbose(stdout);
				fprintf(stdout, "]\n\n");
			}
		} else if (treeSet[0]->numTips != numTaxa) {
			fprintf(stderr, "All trees must have the same number of tips.\n");
			exit(4);
		}
		
		if (maxPartitions == 1) {
			if (partitionLengths[0] != -1) {
				fprintf(stderr, "\nWARNING: The treefile contained partion lengths but only one partition\n");
				fprintf(stderr, "was specified.\n");
			}
			partitionLengths[0] = numSites;
		}

		sumLength = partitionLengths[0];
		i = 1;
		while (sumLength < numSites && i <= maxPartitions) {
			if (!IsTreeAvail(tree_fv)) {
				fprintf(stderr, "\nA set of trees number %d had less partition length (%d) than\n", treeNo + 1, sumLength);
				fprintf(stderr, "was required to make a sequence of length %d.\n", numSites);
				exit(4);
			}
				
			ReadTree(tree_fv, treeSet[i], treeNo+1, treeSet[0]->numTips, treeSet[0]->names, 
						&partitionLengths[i], &partitionRates[i]);
						
			if (treeSet[i]->numTips != numTaxa) {
				fprintf(stderr, "All trees must have the same number of tips.\n");
				exit(4);
			}
			
			sumLength += partitionLengths[i];
			i++;
		}
		if (i > maxPartitions) {
			fprintf(stderr, "\nA set of trees number %d had more partitions (%d) than\n", treeNo + 1, i);
			fprintf(stderr, "was specified in the user options (%d).\n", maxPartitions);
		}
		numPartitions = i;
				
		if (sumLength != numSites) {
			fprintf(stderr, "The sum of the partition lengths in the treefile does not equal\n");
			fprintf(stderr, "the specified number of sites.\n");
			exit(4);
		}
			
		for (i = 0; i < numPartitions; i++)
			CreateSequences(treeSet[i], partitionLengths[i]);
		
		if (numPartitions > 1) {
			sum = 0.0;
			for (i = 0; i < numPartitions; i++)
				sum += partitionRates[i] * partitionLengths[i];
				
			for (i = 0; i < numPartitions; i++)
				partitionRates[i] *= numSites / sum;
		}
		
		if (treeNo==0 && verbose && !quiet) {
			PrintVerbose(stderr);
			InitProgressBar(numTrees*numDatasets);
			DrawProgressBar();
		}

		for (i=0; i<numDatasets; i++) {
			SetCategories();
			
			k = 0;
			for (j = 0; j < numPartitions; j++) {
				scale = partitionRates[j];
				
				if (scaleTrees) { 
					if (!treeSet[j]->rooted) {
						fprintf(stderr, "To scale tree length, they must be rooted and ultrametric.\n");
						exit(4);
					}
					scale *= treeScale/treeSet[j]->totalLength;
				} else if (scaleBranches)
					scale *= branchScale;

				EvolveSequences(treeSet[j], k, partitionLengths[j], scale, ancestor);
				k += partitionLengths[j];
			}
			
			if (writeAncestors)
				WriteAncestralSequences(stdout, treeSet[0]);
			else
				WriteSequences(stdout, (numTrees > 1 ? treeNo+1 : -1), (numDatasets > 1 ? i+1 : -1), treeSet, partitionLengths);

			if (writeRates) {
				WriteRates(stderr);
			}

			if (textFile) {
				while (!feof(text_fv)) {
					ch = fgetc(text_fv);
					if (!feof(text_fv))
						fputc(ch, stdout);
				}
				fputc('\n', stdout);
				rewind(text_fv);
			}
			
			if (verbose && !quiet)
				ProgressBar();
		}
				
		for (i = 0; i < numPartitions; i++)
			DisposeTree(treeSet[i]);
			
		treeNo++;
	} while (IsTreeAvail(tree_fv));
	
/*	for (i = 0; i < maxPartitions; i++)
		FreeTree(treeSet[i]);	*/
	
	if (treeFile)
		fclose(tree_fv);

	if (textFile)
		fclose(text_fv);

	totalSecs = (double)(clock() - totalStart) / CLOCKS_PER_SEC;
	if (!quiet) {
		fprintf(stderr, "Time taken: %G seconds\n", totalSecs);
		if (verboseMemory)
			fprintf(stderr, "Total memory used: %ld\n", totalMem);
	}
	
	return 0;
}
