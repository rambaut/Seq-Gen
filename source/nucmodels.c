/*  
   Sequence Generator - seq-gen, version 1.3.3
   Copyright (c)1996-2011, Andrew Rambaut & Nick Grassly
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
#include <time.h>
#include <math.h>
#include "eigen.h"
#include "nucmodels.h"

char *nucleotides="ACGT";

int nucModel = NUC_NONE;
int equalTstv;

double nucFreq[NUM_NUC] = { 0.25, 0.25, 0.25, 0.25 };
double nucAddFreq[NUM_NUC];
double nucMatrix[MAX_RATE_CATS][SQNUM_NUC];
double nucVector[NUM_NUC];


double freqR, freqY, freqAG, freqCT;
double freqA, freqC, freqG, freqT;
double tstv, kappa;

static double beta, beta_A_R, beta_A_Y;
static double tab1A, tab2A, tab3A;
static double tab1C, tab2C, tab3C;
static double tab1G, tab2G, tab3G;
static double tab1T, tab2T, tab3T;
static double mu, mu_kappa_1;

double nucRelativeRates[NUM_NUC_REL_RATES] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
static double Qij[SQNUM_NUC], Cijk[CUNUM_NUC], Root[NUM_NUC];

void SetHKYMatrix(double *matrix, double len);
void SetHKYVector(double *vector, short state, double len);
void SetF84Matrix(double *matrix, double len);
void SetF84Vector(double *vector, short base, double len);
void SetGTRMatrix(double *matrix, double len);
void SetGTRVector(double *vector, short state, double len);

void CommonMatrix(double aa, double bbR, double bbY, double *matrix);
void CommonVector(double aa, double bbR, double bbY, double *vector, short state);
void CumulativeRows(double *matrix);

void SetupGTR();
void CheckNucFrequencies();

/*************************************/
void SetNucModel(int theNucModel)
{	
	double xi, xv, F84temp1, F84temp2;
	double freqAG, freqCT, freqA2, freqC2, freqG2, freqT2;
	
	nucModel = theNucModel;
	
	freqA=nucFreq[A];
	freqC=nucFreq[C];
	freqG=nucFreq[G];
	freqT=nucFreq[T];
	
	if (nucModel==NUC_F84 || nucModel==NUC_HKY) {
		freqR=freqA+freqG;
		freqY=freqC+freqT;
		freqAG=freqA*freqG;
		freqCT=freqC*freqT;
		
		tab1A=freqA*((1/freqR)-1);
		tab2A=(freqR-freqA)/freqR;
		tab3A=freqA/freqR;
		tab1C=freqC*((1/freqY)-1);
		tab2C=(freqY-freqC)/freqY;
		tab3C=freqC/freqY;
		tab1G=freqG*((1/freqR)-1);
		tab2G=(freqR-freqG)/freqR;
		tab3G=freqG/freqR;
		tab1T=freqT*((1/freqY)-1);
		tab2T=(freqY-freqT)/freqY;
		tab3T=freqT/freqY;
	}
		
	switch (nucModel) {
		case NUC_HKY:
			if (!equalTstv) {
				kappa=(tstv*freqR*freqY)/(freqAG+freqCT);
			} else {
				kappa = 1.0;
				tstv=(kappa*(freqA*freqG + freqC*freqT))/(freqR*freqY);
			}
			beta=-1.0/(2*(freqR*freqY + kappa*(freqAG+freqCT)));
			
			beta_A_R=beta*(1.0+freqR*(kappa-1));
			beta_A_Y=beta*(1.0+freqY*(kappa-1));
		break;
		case NUC_F84:
			freqA2=freqA*freqA;
			freqC2=freqC*freqC;
			freqG2=freqG*freqG;
			freqT2=freqT*freqT;
			F84temp1=freqA2+freqC2+freqG2+freqT2;
			F84temp2=((freqA2/freqR)+(freqC2/freqY)+(freqG2/freqR)+(freqT2/freqY));

			if (!equalTstv) {
				xi=freqR*freqY*(freqR*freqY*tstv-freqAG-freqCT);	
				xv=(freqCT*freqR)+(freqAG*freqY);
				kappa=xi/xv;
			} else {
				kappa = 0.0;
				tstv=(freqY*(freqAG*freqR+freqCT*freqR))/(freqR*freqR*freqY*freqY);
			}	
			mu=-1.0/((1-F84temp1)+(kappa*(1-F84temp2)));
			mu_kappa_1=mu*(kappa+1);
		break;
		case NUC_GTR:
			SetupGTR();
		break;
	}
}

void SetNucMatrix(double *matrix, double len)
{
	switch (nucModel) {
		case NUC_HKY:
			SetHKYMatrix(matrix, len);
		break;
		case NUC_F84:
			SetF84Matrix(matrix, len);
		break;
		case NUC_GTR:
			SetGTRMatrix(matrix, len);
		break;
	}
}

void SetNucVector(double *vector, short state, double len)
{
	switch (nucModel) {
		case NUC_HKY:
			SetHKYVector(vector, state, len);
		break;
		case NUC_F84:
			SetF84Vector(vector, state, len);
		break;
		case NUC_GTR:
			SetGTRVector(vector, state, len);
		break;
	}
}

void SetHKYMatrix(double *matrix, double len)
{	
	double aa, bbR, bbY;
	
	aa=exp(beta*len);
	bbR=exp(beta_A_R*len);
	bbY=exp(beta_A_Y*len);
	
	CommonMatrix(aa, bbR, bbY, matrix);
}

void SetHKYVector(double *vector, short state, double len)
{	
	double aa, bbR, bbY;

	aa=exp(beta*len);
	bbR=exp(beta_A_R*len);
	bbY=exp(beta_A_Y*len);

	CommonVector(aa, bbR, bbY, vector, state);
}

void SetF84Matrix(double *matrix, double len)
{	
	double aa, bbR, bbY;
	
	aa=exp(mu*len);
	bbR=bbY=exp(mu_kappa_1*len);

	CommonMatrix(aa, bbR, bbY, matrix);
}

void SetF84Vector(double *vector, short state, double len)
{	
	double aa, bbR, bbY;

	aa=exp(mu*len);
	bbR=bbY=exp(mu_kappa_1*len);
	
	CommonVector(aa, bbR, bbY, vector, state);
}

void SetGTRMatrix(double *matrix, double len)
{	
	int i,j,k;
	double expt[4];
	double *P;
	
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
	P=matrix;
	if (len<1e-6) { 
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				if (i==j)
					*P=1.0;
				else 	
					*P=0.0;
				P++;
			}
		}
		return; 
	}
	
	for (k=1; k<4; k++) {
		expt[k]=exp(len*Root[k]);
	}
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			(*P)=Cijk[i*4*4+j*4+0];
			for (k=1; k<4; k++) {
				(*P)+=Cijk[i*4*4+j*4+k]*expt[k];
			}
			P++;
		}
	}
	
	CumulativeRows(matrix);
}

void SetGTRVector(double *vector, short state, double len)
{	
	int i,j,k;
	double expt[4];
	double *P;

	P=vector;
	if (len<1e-6) { 
		for (i=0; i<4; i++) {
			if (i==state)
				*P=1.0;
			else 	
				*P=0.0;
			P++;
		}
		return; 
	}
	for (k=1; k<4; k++) {
		expt[k]=exp(len*Root[k]);
	}
	
	for (j=0; j<4; j++) {
		(*P)=Cijk[state*4*4+j*4+0];
		for (k=1; k<4; k++) {
			(*P)+=Cijk[state*4*4+j*4+k]*expt[k];
		}
		P++;
	}
	
	vector[1]+=vector[0];
	vector[2]+=vector[1];
	vector[3]+=vector[2];
}

		
#define PIJ_SAME_A freqA+(tab1A*aa)+(tab2A*bbR)
#define PIJ_TS_A freqA+(tab1A*aa)-(tab3A*bbR)
#define PIJ_TV_A freqA*(1-aa)

#define PIJ_SAME_C freqC+(tab1C*aa)+(tab2C*bbY)
#define PIJ_TS_C freqC+(tab1C*aa)-(tab3C*bbY)
#define PIJ_TV_C freqC*(1-aa)
	
#define PIJ_SAME_G freqG+(tab1G*aa)+(tab2G*bbR)
#define PIJ_TS_G freqG+(tab1G*aa)-(tab3G*bbR)
#define PIJ_TV_G freqG*(1-aa)
	
#define PIJ_SAME_T freqT+(tab1T*aa)+(tab2T*bbY)
#define PIJ_TS_T freqT+(tab1T*aa)-(tab3T*bbY)
#define PIJ_TV_T freqT*(1-aa)	
	
void CommonMatrix(double aa, double bbR, double bbY, double *matrix) 
{
	matrix[0]=PIJ_SAME_A;
	matrix[1]=PIJ_TV_C;
	matrix[2]=PIJ_TS_G;
	matrix[3]=PIJ_TV_T;

	matrix[4]=PIJ_TV_A;
	matrix[5]=PIJ_SAME_C;
	matrix[6]=PIJ_TV_G;
	matrix[7]=PIJ_TS_T;
	
	matrix[8]=PIJ_TS_A;
	matrix[9]=matrix[1];  /* PIJ_TV_C */
	matrix[10]=PIJ_SAME_G;
	matrix[11]=matrix[3]; /* PIJ_TV_T */
	
	matrix[12]=matrix[4]; /* PIJ_TV_A */
	matrix[13]=PIJ_TS_C;
	matrix[14]=matrix[6]; /* PIJ_TV_G */
	matrix[15]=PIJ_SAME_T;

	CumulativeRows(matrix);
}

void CumulativeRows(double *matrix) 
{
/* the rows are cumulative to help with picking one using
   a random number */
	matrix[1]+=matrix[0];
	matrix[2]+=matrix[1];
	matrix[3]+=matrix[2]; /* This should always be 1.0... */

	matrix[5]+=matrix[4];
	matrix[6]+=matrix[5];
	matrix[7]+=matrix[6]; /* ...but it is easier to spot bugs... */
	
	matrix[9]+=matrix[8];
	matrix[10]+=matrix[9];
	matrix[11]+=matrix[10]; /* ...though less efficient... */
	
	matrix[13]+=matrix[12];
	matrix[14]+=matrix[13];
	matrix[15]+=matrix[14]; /* ...but probably not much. */
}


void CommonVector(double aa, double bbR, double bbY, double *vector, short state)
{	
	switch (state) {
		case 0:
			vector[0]=PIJ_SAME_A;
			vector[1]=PIJ_TV_C+vector[0];
			vector[2]=PIJ_TS_G+vector[1];
			vector[3]=PIJ_TV_T+vector[2];
		break;
		case 1:
			vector[0]=PIJ_TV_A;
			vector[1]=PIJ_SAME_C+vector[0];
			vector[2]=PIJ_TV_G+vector[1];
			vector[3]=PIJ_TS_T+vector[2];
		break;
		case 2:
			vector[0]=PIJ_TS_A;
			vector[1]=PIJ_TV_C+vector[0];
			vector[2]=PIJ_SAME_G+vector[1];
			vector[3]=PIJ_TV_T+vector[2];
		break;
		case 3:
			vector[0]=PIJ_TV_A;
			vector[1]=PIJ_TS_C+vector[0];
			vector[2]=PIJ_TV_G+vector[1];
			vector[3]=PIJ_SAME_T+vector[2];
		break;
	}
}

void SetupGTR()
{
	int i,j,k;
	double mr;
	double sum;
	double U[SQNUM_NUC], V[SQNUM_NUC], T1[SQNUM_NUC], T2[SQNUM_NUC];

	CheckNucFrequencies();
	
	k=0;
	for (i=0; i<NUM_NUC-1; i++) {
		for (j=i+1; j<NUM_NUC; j++) {
			Qij[i*NUM_NUC+j] = Qij[j*NUM_NUC+i] = nucRelativeRates[k++];
		}
	}
	
	for (i=0; i<NUM_NUC; i++) {
		for (j=0; j<NUM_NUC; j++) { 
			Qij[i*NUM_NUC+j] *= nucFreq[j];
		}
	}
		
	mr=0;		
	for (i=0; i<NUM_NUC; i++) {
		sum = 0;
		Qij[i*NUM_NUC+i]=0; 
		for (j=0; j<NUM_NUC; j++) { 
			sum += Qij[i*NUM_NUC+j];
		}
		Qij[i*NUM_NUC+i] = -sum; 
		mr += nucFreq[i] * sum;
	}
	
	abyx(1.0/mr, Qij, SQNUM_NUC);

	if ((k=eigen(1, Qij, NUM_NUC, Root, T1, U, V, T2))!=0) {
		fprintf(stderr, "\ncomplex roots in SetupGTR");
		exit(0);
	}
	xtoy (U, V, SQNUM_NUC);
	matinv (V, NUM_NUC, NUM_NUC, T1);
	for (i=0; i<NUM_NUC; i++) {
   		for (j=0; j<NUM_NUC; j++) {
   			for (k=0; k<NUM_NUC; k++) {
   				Cijk[i*SQNUM_NUC+j*NUM_NUC+k] = U[i*NUM_NUC+k]*V[k*NUM_NUC+j];
   			}
   		}
   	}
}

/**
 * Ensures that frequencies are not smaller than MINFREQ and
 * that two frequencies differ by at least 2*MINFDIFF.
 * This avoids potential problems later when eigenvalues
 * are computed.
 */
void CheckNucFrequencies()
{
	int i, j;
	double diff;
	
	// required frequency difference
	double MINFDIFF = 1.0E-10;

	// lower limit on frequency
	double MINFREQ = 1.0E-10;

	int maxi = 0;
	double sum = 0.0;
	double maxfreq = 0.0;
	for (i = 0; i < NUM_NUC; i++) {
		if (nucFreq[i] < MINFREQ) nucFreq[i] = MINFREQ;
		if (nucFreq[i] > maxfreq) {
			maxfreq = nucFreq[i];
			maxi = i;
		}
		sum += nucFreq[i];
	}
	
	diff = 1.0 - sum;
	nucFreq[maxi] += diff;

	for (i = 0; i < NUM_NUC - 1; i++) {
		for (j = i+1; j < NUM_NUC; j++) {
			if (nucFreq[i] == nucFreq[j]) {
				nucFreq[i] += MINFDIFF;
				nucFreq[j] -= MINFDIFF;
			}
		}
	}
}

