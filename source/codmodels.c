/*
 *  codmodels.c
 *  
 *
 *  Created by Daniel Wilson on 3/6/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "paml.h"
#include "codmodels.h"
#include "paml.h"

/* # Some R code to generate the codons string
base = c("T","C","A","G")
l=1; str = ""
for(i in 1:4) {
  for(j in 1:4) {
    for(k in 1:4) {
      if(l!=11 & l!=12 & l!=15) {
        str = paste(str,base[i],base[j],base[k],sep="")
      }
      l = l+1;
    }
  }
}*/
char* codons="TTTTTCTTATTGTCTTCCTCATCGTATTACTGTTGCTGGCTTCTCCTACTGCCTCCCCCACCGCATCACCAACAGCGTCGCCGACGGATTATCATAATGACTACCACAACGAATAACAAAAAGAGTAGCAGAAGGGTTGTCGTAGTGGCTGCCGCAGCGGATGACGAAGAGGGTGGCGGAGGG";

int codFreqSet;

int codModel = COD_NONE;
double codFreq[NUM_COD];
double sqrtCodFreq[NUM_COD];
double codAddFreq[NUM_COD];
double codMatrix[MAX_RATE_CATS][SQNUM_COD];
double codVector[NUM_COD];

double tstv, kappa,omega;

static double Qij[SQNUM_COD], Cijk[CUNUM_COD], Root[NUM_COD];

void SetupCodMatrix();
void SetupNY98Matrix();
void SetCodFrequencies(double* inFrequencies);
void CheckCodFrequencies();

static double equalFrequencies[NUM_COD] = {0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508, 0.0163934426229508};

/*************************************/
void SetCodModel(int theCodModel) {
	codModel = theCodModel;

	kappa = tstv;

	if(!codFreqSet) {
		SetCodFrequencies(equalFrequencies);
	}
	else {
		CheckCodFrequencies();
	}
	
	SetupCodMatrix();
}

void SetCodFrequencies(double *inFrequencies)
{
	int i;
	for (i=0; i<NUM_COD; i++) {
		codFreq[i] = inFrequencies[i];
		sqrtCodFreq[i] = sqrt(codFreq[i]);
	}
}


void SetupCodMatrix() {
	int i,j,k;
	double mr;
	double sum;
	double U[SQNUM_COD], T1[SQNUM_COD];

	switch(codModel) {
		case COD_NY98:
			SetupNY98Matrix();
		break;
		default:
			fprintf(stderr,"\nModel not specified in SetupCodMatrix()");
			exit(0);
	}

	/* Multiply by stationary frequencies - no longer same routine
	for (i=0; i<NUM_COD; i++) {
		for (j=0; j<NUM_COD; j++) { 
			Qij[i*NUM_COD+j] *= codFreq[j];
		}
	}*/

	/* Compute the diagonal and mean rate */
	mr=0;		
	for (i=0; i<NUM_COD; i++) {
		sum = 0;
		Qij[i*NUM_COD+i]=0; 
		for (j=0; j<NUM_COD; j++) { 
			sum += Qij[i*NUM_COD+j];
			/* New method to calculate mean rate */
			if(j>i) mr += sqrtCodFreq[i]*sqrtCodFreq[j]*Qij[i*NUM_COD+j];
		}
		Qij[i*NUM_COD+i] = -sum; 
		/*mr += codFreq[i] * sum;*/
	}
	mr *= 2.0;

	/* Divide through by the average rate */
	abyx(1.0/mr, Qij, SQNUM_COD);

	/* Output rate matrix for inspection
	int ii,jj;
	FILE *Qijfile = fopen("Qij.txt","w");
	for(jj=0;jj<61;jj++) {
		for(ii=0;ii<61;ii++) {
			fprintf(Qijfile,"%g ",Qij[ii*NUM_COD+jj]);
		}
	}
	fclose(Qijfile);*/

	/* Replace by Real Symmetric Algorithm - more stable */
	xtoy(Qij,U,SQNUM_COD);
	eigenRealSym(U,NUM_COD,Root,T1);
	/* This routines was failing without flagging an error
	if ((k=eigen(1, Qij, NUM_COD, Root, T1, U, V, T2))!=0) {
		fprintf(stderr, "\ncomplex roots in SetupCodMatrix");
		exit(0);
	}*/
	/* U is symmetric by definition. matinv was also failing with flagging an error
	xtoy (U, V, SQNUM_COD);
	matinv (V, NUM_COD, NUM_COD, T1);*/

	/* Output eigenvalues and vectors for inspection
	FILE *eigUfile = fopen("eigUfile.txt","w");
	FILE *eigVfile = fopen("eigVfile.txt","w");
	FILE *eigDfile = fopen("eigDfile.txt","w");
	for(jj=0;jj<61;jj++) {
		for(ii=0;ii<61;ii++) {
			fprintf(eigUfile,"%g ",U[ii*NUM_COD+jj]);
			fprintf(eigVfile,"%g ",U[jj*NUM_COD+ii]);
		}
		fprintf(eigDfile,"%g ",Root[jj]);
	}
	fclose(eigUfile);
	fclose(eigVfile);
	fclose(eigDfile);*/
	
	/*Check matrix inversion
	for(i=0; i<NUM_COD; i++) {
		for(j=0; j<NUM_COD; j++) {
			double el = 0;
			for(k=0; k<NUM_COD; k++) {
				el += U[i*NUM_COD+k] * U[j*NUM_COD+k];
			}
			if((i==j & fabs(el-1.)>1e-6) ||
			   (i!=j & fabs(el)>1e-6)) {
				fprintf(stderr,"\nI[%d][%d]=%g\n\n",i,j,el);
				exit(0);
			}
		}
	}*/
	
	for (i=0; i<NUM_COD; i++) {
   		for (j=0; j<NUM_COD; j++) {
   			for (k=0; k<NUM_COD; k++) {
				/* New method for calculating Cijk */
   				Cijk[i*SQNUM_COD+j*NUM_COD+k] = U[i*NUM_COD+k]*U[j*NUM_COD+k]*sqrtCodFreq[j]/sqrtCodFreq[i];
   			}
   		}
   	}
	
}

void SetupNY98Matrix() {
	int i,j;
	/* Initialize to zero */
	for(i=0;i<NUM_COD;i++) for(j=0;j<NUM_COD;j++) Qij[i*NUM_COD+j] = 0;

	/* Fill in upper triangle */
	Qij[0*NUM_COD+1] = sqrtCodFreq[0]*sqrtCodFreq[1]*kappa;
	Qij[0*NUM_COD+2] = sqrtCodFreq[0]*sqrtCodFreq[2]*omega;
	Qij[0*NUM_COD+3] = sqrtCodFreq[0]*sqrtCodFreq[3]*omega;
	Qij[0*NUM_COD+4] = sqrtCodFreq[0]*sqrtCodFreq[4]*kappa*omega;
	Qij[0*NUM_COD+8] = sqrtCodFreq[0]*sqrtCodFreq[8]*omega;
	Qij[0*NUM_COD+10] = sqrtCodFreq[0]*sqrtCodFreq[10]*omega;
	Qij[0*NUM_COD+13] = sqrtCodFreq[0]*sqrtCodFreq[13]*kappa*omega;
	Qij[0*NUM_COD+29] = sqrtCodFreq[0]*sqrtCodFreq[29]*omega;
	Qij[0*NUM_COD+45] = sqrtCodFreq[0]*sqrtCodFreq[45]*omega;
	Qij[1*NUM_COD+2] = sqrtCodFreq[1]*sqrtCodFreq[2]*omega;
	Qij[1*NUM_COD+3] = sqrtCodFreq[1]*sqrtCodFreq[3]*omega;
	Qij[1*NUM_COD+5] = sqrtCodFreq[1]*sqrtCodFreq[5]*kappa*omega;
	Qij[1*NUM_COD+9] = sqrtCodFreq[1]*sqrtCodFreq[9]*omega;
	Qij[1*NUM_COD+11] = sqrtCodFreq[1]*sqrtCodFreq[11]*omega;
	Qij[1*NUM_COD+14] = sqrtCodFreq[1]*sqrtCodFreq[14]*kappa*omega;
	Qij[1*NUM_COD+30] = sqrtCodFreq[1]*sqrtCodFreq[30]*omega;
	Qij[1*NUM_COD+46] = sqrtCodFreq[1]*sqrtCodFreq[46]*omega;
	Qij[2*NUM_COD+3] = sqrtCodFreq[2]*sqrtCodFreq[3]*kappa;
	Qij[2*NUM_COD+6] = sqrtCodFreq[2]*sqrtCodFreq[6]*kappa*omega;
	Qij[2*NUM_COD+15] = sqrtCodFreq[2]*sqrtCodFreq[15]*kappa;
	Qij[2*NUM_COD+31] = sqrtCodFreq[2]*sqrtCodFreq[31]*omega;
	Qij[2*NUM_COD+47] = sqrtCodFreq[2]*sqrtCodFreq[47]*omega;
	Qij[3*NUM_COD+7] = sqrtCodFreq[3]*sqrtCodFreq[7]*kappa*omega;
	Qij[3*NUM_COD+12] = sqrtCodFreq[3]*sqrtCodFreq[12]*omega;
	Qij[3*NUM_COD+16] = sqrtCodFreq[3]*sqrtCodFreq[16]*kappa;
	Qij[3*NUM_COD+32] = sqrtCodFreq[3]*sqrtCodFreq[32]*omega;
	Qij[3*NUM_COD+48] = sqrtCodFreq[3]*sqrtCodFreq[48]*omega;
	Qij[4*NUM_COD+5] = sqrtCodFreq[4]*sqrtCodFreq[5]*kappa;
	Qij[4*NUM_COD+6] = sqrtCodFreq[4]*sqrtCodFreq[6];
	Qij[4*NUM_COD+7] = sqrtCodFreq[4]*sqrtCodFreq[7];
	Qij[4*NUM_COD+8] = sqrtCodFreq[4]*sqrtCodFreq[8]*omega;
	Qij[4*NUM_COD+10] = sqrtCodFreq[4]*sqrtCodFreq[10]*omega;
	Qij[4*NUM_COD+17] = sqrtCodFreq[4]*sqrtCodFreq[17]*kappa*omega;
	Qij[4*NUM_COD+33] = sqrtCodFreq[4]*sqrtCodFreq[33]*omega;
	Qij[4*NUM_COD+49] = sqrtCodFreq[4]*sqrtCodFreq[49]*omega;
	Qij[5*NUM_COD+6] = sqrtCodFreq[5]*sqrtCodFreq[6];
	Qij[5*NUM_COD+7] = sqrtCodFreq[5]*sqrtCodFreq[7];
	Qij[5*NUM_COD+9] = sqrtCodFreq[5]*sqrtCodFreq[9]*omega;
	Qij[5*NUM_COD+11] = sqrtCodFreq[5]*sqrtCodFreq[11]*omega;
	Qij[5*NUM_COD+18] = sqrtCodFreq[5]*sqrtCodFreq[18]*kappa*omega;
	Qij[5*NUM_COD+34] = sqrtCodFreq[5]*sqrtCodFreq[34]*omega;
	Qij[5*NUM_COD+50] = sqrtCodFreq[5]*sqrtCodFreq[50]*omega;
	Qij[6*NUM_COD+7] = sqrtCodFreq[6]*sqrtCodFreq[7]*kappa;
	Qij[6*NUM_COD+19] = sqrtCodFreq[6]*sqrtCodFreq[19]*kappa*omega;
	Qij[6*NUM_COD+35] = sqrtCodFreq[6]*sqrtCodFreq[35]*omega;
	Qij[6*NUM_COD+51] = sqrtCodFreq[6]*sqrtCodFreq[51]*omega;
	Qij[7*NUM_COD+12] = sqrtCodFreq[7]*sqrtCodFreq[12]*omega;
	Qij[7*NUM_COD+20] = sqrtCodFreq[7]*sqrtCodFreq[20]*kappa*omega;
	Qij[7*NUM_COD+36] = sqrtCodFreq[7]*sqrtCodFreq[36]*omega;
	Qij[7*NUM_COD+52] = sqrtCodFreq[7]*sqrtCodFreq[52]*omega;
	Qij[8*NUM_COD+9] = sqrtCodFreq[8]*sqrtCodFreq[9]*kappa;
	Qij[8*NUM_COD+10] = sqrtCodFreq[8]*sqrtCodFreq[10]*kappa*omega;
	Qij[8*NUM_COD+21] = sqrtCodFreq[8]*sqrtCodFreq[21]*kappa*omega;
	Qij[8*NUM_COD+37] = sqrtCodFreq[8]*sqrtCodFreq[37]*omega;
	Qij[8*NUM_COD+53] = sqrtCodFreq[8]*sqrtCodFreq[53]*omega;
	Qij[9*NUM_COD+11] = sqrtCodFreq[9]*sqrtCodFreq[11]*kappa*omega;
	Qij[9*NUM_COD+22] = sqrtCodFreq[9]*sqrtCodFreq[22]*kappa*omega;
	Qij[9*NUM_COD+38] = sqrtCodFreq[9]*sqrtCodFreq[38]*omega;
	Qij[9*NUM_COD+54] = sqrtCodFreq[9]*sqrtCodFreq[54]*omega;
	Qij[10*NUM_COD+11] = sqrtCodFreq[10]*sqrtCodFreq[11]*kappa;
	Qij[10*NUM_COD+12] = sqrtCodFreq[10]*sqrtCodFreq[12]*omega;
	Qij[10*NUM_COD+25] = sqrtCodFreq[10]*sqrtCodFreq[25]*kappa*omega;
	Qij[10*NUM_COD+41] = sqrtCodFreq[10]*sqrtCodFreq[41]*omega;
	Qij[10*NUM_COD+57] = sqrtCodFreq[10]*sqrtCodFreq[57]*omega;
	Qij[11*NUM_COD+12] = sqrtCodFreq[11]*sqrtCodFreq[12]*omega;
	Qij[11*NUM_COD+26] = sqrtCodFreq[11]*sqrtCodFreq[26]*kappa*omega;
	Qij[11*NUM_COD+42] = sqrtCodFreq[11]*sqrtCodFreq[42]*omega;
	Qij[11*NUM_COD+58] = sqrtCodFreq[11]*sqrtCodFreq[58]*omega;
	Qij[12*NUM_COD+28] = sqrtCodFreq[12]*sqrtCodFreq[28]*kappa*omega;
	Qij[12*NUM_COD+44] = sqrtCodFreq[12]*sqrtCodFreq[44]*omega;
	Qij[12*NUM_COD+60] = sqrtCodFreq[12]*sqrtCodFreq[60]*omega;
	Qij[13*NUM_COD+14] = sqrtCodFreq[13]*sqrtCodFreq[14]*kappa;
	Qij[13*NUM_COD+15] = sqrtCodFreq[13]*sqrtCodFreq[15];
	Qij[13*NUM_COD+16] = sqrtCodFreq[13]*sqrtCodFreq[16];
	Qij[13*NUM_COD+17] = sqrtCodFreq[13]*sqrtCodFreq[17]*kappa*omega;
	Qij[13*NUM_COD+21] = sqrtCodFreq[13]*sqrtCodFreq[21]*omega;
	Qij[13*NUM_COD+25] = sqrtCodFreq[13]*sqrtCodFreq[25]*omega;
	Qij[13*NUM_COD+29] = sqrtCodFreq[13]*sqrtCodFreq[29]*omega;
	Qij[13*NUM_COD+45] = sqrtCodFreq[13]*sqrtCodFreq[45]*omega;
	Qij[14*NUM_COD+15] = sqrtCodFreq[14]*sqrtCodFreq[15];
	Qij[14*NUM_COD+16] = sqrtCodFreq[14]*sqrtCodFreq[16];
	Qij[14*NUM_COD+18] = sqrtCodFreq[14]*sqrtCodFreq[18]*kappa*omega;
	Qij[14*NUM_COD+22] = sqrtCodFreq[14]*sqrtCodFreq[22]*omega;
	Qij[14*NUM_COD+26] = sqrtCodFreq[14]*sqrtCodFreq[26]*omega;
	Qij[14*NUM_COD+30] = sqrtCodFreq[14]*sqrtCodFreq[30]*omega;
	Qij[14*NUM_COD+46] = sqrtCodFreq[14]*sqrtCodFreq[46]*omega;
	Qij[15*NUM_COD+16] = sqrtCodFreq[15]*sqrtCodFreq[16]*kappa;
	Qij[15*NUM_COD+19] = sqrtCodFreq[15]*sqrtCodFreq[19]*kappa*omega;
	Qij[15*NUM_COD+23] = sqrtCodFreq[15]*sqrtCodFreq[23]*omega;
	Qij[15*NUM_COD+27] = sqrtCodFreq[15]*sqrtCodFreq[27]*omega;
	Qij[15*NUM_COD+31] = sqrtCodFreq[15]*sqrtCodFreq[31]*omega;
	Qij[15*NUM_COD+47] = sqrtCodFreq[15]*sqrtCodFreq[47]*omega;
	Qij[16*NUM_COD+20] = sqrtCodFreq[16]*sqrtCodFreq[20]*kappa*omega;
	Qij[16*NUM_COD+24] = sqrtCodFreq[16]*sqrtCodFreq[24]*omega;
	Qij[16*NUM_COD+28] = sqrtCodFreq[16]*sqrtCodFreq[28]*omega;
	Qij[16*NUM_COD+32] = sqrtCodFreq[16]*sqrtCodFreq[32]*omega;
	Qij[16*NUM_COD+48] = sqrtCodFreq[16]*sqrtCodFreq[48]*omega;
	Qij[17*NUM_COD+18] = sqrtCodFreq[17]*sqrtCodFreq[18]*kappa;
	Qij[17*NUM_COD+19] = sqrtCodFreq[17]*sqrtCodFreq[19];
	Qij[17*NUM_COD+20] = sqrtCodFreq[17]*sqrtCodFreq[20];
	Qij[17*NUM_COD+21] = sqrtCodFreq[17]*sqrtCodFreq[21]*omega;
	Qij[17*NUM_COD+25] = sqrtCodFreq[17]*sqrtCodFreq[25]*omega;
	Qij[17*NUM_COD+33] = sqrtCodFreq[17]*sqrtCodFreq[33]*omega;
	Qij[17*NUM_COD+49] = sqrtCodFreq[17]*sqrtCodFreq[49]*omega;
	Qij[18*NUM_COD+19] = sqrtCodFreq[18]*sqrtCodFreq[19];
	Qij[18*NUM_COD+20] = sqrtCodFreq[18]*sqrtCodFreq[20];
	Qij[18*NUM_COD+22] = sqrtCodFreq[18]*sqrtCodFreq[22]*omega;
	Qij[18*NUM_COD+26] = sqrtCodFreq[18]*sqrtCodFreq[26]*omega;
	Qij[18*NUM_COD+34] = sqrtCodFreq[18]*sqrtCodFreq[34]*omega;
	Qij[18*NUM_COD+50] = sqrtCodFreq[18]*sqrtCodFreq[50]*omega;
	Qij[19*NUM_COD+20] = sqrtCodFreq[19]*sqrtCodFreq[20]*kappa;
	Qij[19*NUM_COD+23] = sqrtCodFreq[19]*sqrtCodFreq[23]*omega;
	Qij[19*NUM_COD+27] = sqrtCodFreq[19]*sqrtCodFreq[27]*omega;
	Qij[19*NUM_COD+35] = sqrtCodFreq[19]*sqrtCodFreq[35]*omega;
	Qij[19*NUM_COD+51] = sqrtCodFreq[19]*sqrtCodFreq[51]*omega;
	Qij[20*NUM_COD+24] = sqrtCodFreq[20]*sqrtCodFreq[24]*omega;
	Qij[20*NUM_COD+28] = sqrtCodFreq[20]*sqrtCodFreq[28]*omega;
	Qij[20*NUM_COD+36] = sqrtCodFreq[20]*sqrtCodFreq[36]*omega;
	Qij[20*NUM_COD+52] = sqrtCodFreq[20]*sqrtCodFreq[52]*omega;
	Qij[21*NUM_COD+22] = sqrtCodFreq[21]*sqrtCodFreq[22]*kappa;
	Qij[21*NUM_COD+23] = sqrtCodFreq[21]*sqrtCodFreq[23]*omega;
	Qij[21*NUM_COD+24] = sqrtCodFreq[21]*sqrtCodFreq[24]*omega;
	Qij[21*NUM_COD+25] = sqrtCodFreq[21]*sqrtCodFreq[25]*kappa*omega;
	Qij[21*NUM_COD+37] = sqrtCodFreq[21]*sqrtCodFreq[37]*omega;
	Qij[21*NUM_COD+53] = sqrtCodFreq[21]*sqrtCodFreq[53]*omega;
	Qij[22*NUM_COD+23] = sqrtCodFreq[22]*sqrtCodFreq[23]*omega;
	Qij[22*NUM_COD+24] = sqrtCodFreq[22]*sqrtCodFreq[24]*omega;
	Qij[22*NUM_COD+26] = sqrtCodFreq[22]*sqrtCodFreq[26]*kappa*omega;
	Qij[22*NUM_COD+38] = sqrtCodFreq[22]*sqrtCodFreq[38]*omega;
	Qij[22*NUM_COD+54] = sqrtCodFreq[22]*sqrtCodFreq[54]*omega;
	Qij[23*NUM_COD+24] = sqrtCodFreq[23]*sqrtCodFreq[24]*kappa;
	Qij[23*NUM_COD+27] = sqrtCodFreq[23]*sqrtCodFreq[27]*kappa*omega;
	Qij[23*NUM_COD+39] = sqrtCodFreq[23]*sqrtCodFreq[39]*omega;
	Qij[23*NUM_COD+55] = sqrtCodFreq[23]*sqrtCodFreq[55]*omega;
	Qij[24*NUM_COD+28] = sqrtCodFreq[24]*sqrtCodFreq[28]*kappa*omega;
	Qij[24*NUM_COD+40] = sqrtCodFreq[24]*sqrtCodFreq[40]*omega;
	Qij[24*NUM_COD+56] = sqrtCodFreq[24]*sqrtCodFreq[56]*omega;
	Qij[25*NUM_COD+26] = sqrtCodFreq[25]*sqrtCodFreq[26]*kappa;
	Qij[25*NUM_COD+27] = sqrtCodFreq[25]*sqrtCodFreq[27];
	Qij[25*NUM_COD+28] = sqrtCodFreq[25]*sqrtCodFreq[28];
	Qij[25*NUM_COD+41] = sqrtCodFreq[25]*sqrtCodFreq[41]*omega;
	Qij[25*NUM_COD+57] = sqrtCodFreq[25]*sqrtCodFreq[57]*omega;
	Qij[26*NUM_COD+27] = sqrtCodFreq[26]*sqrtCodFreq[27];
	Qij[26*NUM_COD+28] = sqrtCodFreq[26]*sqrtCodFreq[28];
	Qij[26*NUM_COD+42] = sqrtCodFreq[26]*sqrtCodFreq[42]*omega;
	Qij[26*NUM_COD+58] = sqrtCodFreq[26]*sqrtCodFreq[58]*omega;
	Qij[27*NUM_COD+28] = sqrtCodFreq[27]*sqrtCodFreq[28]*kappa;
	Qij[27*NUM_COD+43] = sqrtCodFreq[27]*sqrtCodFreq[43];
	Qij[27*NUM_COD+59] = sqrtCodFreq[27]*sqrtCodFreq[59]*omega;
	Qij[28*NUM_COD+44] = sqrtCodFreq[28]*sqrtCodFreq[44];
	Qij[28*NUM_COD+60] = sqrtCodFreq[28]*sqrtCodFreq[60]*omega;
	Qij[29*NUM_COD+30] = sqrtCodFreq[29]*sqrtCodFreq[30]*kappa;
	Qij[29*NUM_COD+31] = sqrtCodFreq[29]*sqrtCodFreq[31];
	Qij[29*NUM_COD+32] = sqrtCodFreq[29]*sqrtCodFreq[32]*omega;
	Qij[29*NUM_COD+33] = sqrtCodFreq[29]*sqrtCodFreq[33]*kappa*omega;
	Qij[29*NUM_COD+37] = sqrtCodFreq[29]*sqrtCodFreq[37]*omega;
	Qij[29*NUM_COD+41] = sqrtCodFreq[29]*sqrtCodFreq[41]*omega;
	Qij[29*NUM_COD+45] = sqrtCodFreq[29]*sqrtCodFreq[45]*kappa*omega;
	Qij[30*NUM_COD+31] = sqrtCodFreq[30]*sqrtCodFreq[31];
	Qij[30*NUM_COD+32] = sqrtCodFreq[30]*sqrtCodFreq[32]*omega;
	Qij[30*NUM_COD+34] = sqrtCodFreq[30]*sqrtCodFreq[34]*kappa*omega;
	Qij[30*NUM_COD+38] = sqrtCodFreq[30]*sqrtCodFreq[38]*omega;
	Qij[30*NUM_COD+42] = sqrtCodFreq[30]*sqrtCodFreq[42]*omega;
	Qij[30*NUM_COD+46] = sqrtCodFreq[30]*sqrtCodFreq[46]*kappa*omega;
	Qij[31*NUM_COD+32] = sqrtCodFreq[31]*sqrtCodFreq[32]*kappa*omega;
	Qij[31*NUM_COD+35] = sqrtCodFreq[31]*sqrtCodFreq[35]*kappa*omega;
	Qij[31*NUM_COD+39] = sqrtCodFreq[31]*sqrtCodFreq[39]*omega;
	Qij[31*NUM_COD+43] = sqrtCodFreq[31]*sqrtCodFreq[43]*omega;
	Qij[31*NUM_COD+47] = sqrtCodFreq[31]*sqrtCodFreq[47]*kappa*omega;
	Qij[32*NUM_COD+36] = sqrtCodFreq[32]*sqrtCodFreq[36]*kappa*omega;
	Qij[32*NUM_COD+40] = sqrtCodFreq[32]*sqrtCodFreq[40]*omega;
	Qij[32*NUM_COD+44] = sqrtCodFreq[32]*sqrtCodFreq[44]*omega;
	Qij[32*NUM_COD+48] = sqrtCodFreq[32]*sqrtCodFreq[48]*kappa*omega;
	Qij[33*NUM_COD+34] = sqrtCodFreq[33]*sqrtCodFreq[34]*kappa;
	Qij[33*NUM_COD+35] = sqrtCodFreq[33]*sqrtCodFreq[35];
	Qij[33*NUM_COD+36] = sqrtCodFreq[33]*sqrtCodFreq[36];
	Qij[33*NUM_COD+37] = sqrtCodFreq[33]*sqrtCodFreq[37]*omega;
	Qij[33*NUM_COD+41] = sqrtCodFreq[33]*sqrtCodFreq[41]*omega;
	Qij[33*NUM_COD+49] = sqrtCodFreq[33]*sqrtCodFreq[49]*kappa*omega;
	Qij[34*NUM_COD+35] = sqrtCodFreq[34]*sqrtCodFreq[35];
	Qij[34*NUM_COD+36] = sqrtCodFreq[34]*sqrtCodFreq[36];
	Qij[34*NUM_COD+38] = sqrtCodFreq[34]*sqrtCodFreq[38]*omega;
	Qij[34*NUM_COD+42] = sqrtCodFreq[34]*sqrtCodFreq[42]*omega;
	Qij[34*NUM_COD+50] = sqrtCodFreq[34]*sqrtCodFreq[50]*kappa*omega;
	Qij[35*NUM_COD+36] = sqrtCodFreq[35]*sqrtCodFreq[36]*kappa;
	Qij[35*NUM_COD+39] = sqrtCodFreq[35]*sqrtCodFreq[39]*omega;
	Qij[35*NUM_COD+43] = sqrtCodFreq[35]*sqrtCodFreq[43]*omega;
	Qij[35*NUM_COD+51] = sqrtCodFreq[35]*sqrtCodFreq[51]*kappa*omega;
	Qij[36*NUM_COD+40] = sqrtCodFreq[36]*sqrtCodFreq[40]*omega;
	Qij[36*NUM_COD+44] = sqrtCodFreq[36]*sqrtCodFreq[44]*omega;
	Qij[36*NUM_COD+52] = sqrtCodFreq[36]*sqrtCodFreq[52]*kappa*omega;
	Qij[37*NUM_COD+38] = sqrtCodFreq[37]*sqrtCodFreq[38]*kappa;
	Qij[37*NUM_COD+39] = sqrtCodFreq[37]*sqrtCodFreq[39]*omega;
	Qij[37*NUM_COD+40] = sqrtCodFreq[37]*sqrtCodFreq[40]*omega;
	Qij[37*NUM_COD+41] = sqrtCodFreq[37]*sqrtCodFreq[41]*kappa*omega;
	Qij[37*NUM_COD+53] = sqrtCodFreq[37]*sqrtCodFreq[53]*kappa*omega;
	Qij[38*NUM_COD+39] = sqrtCodFreq[38]*sqrtCodFreq[39]*omega;
	Qij[38*NUM_COD+40] = sqrtCodFreq[38]*sqrtCodFreq[40]*omega;
	Qij[38*NUM_COD+42] = sqrtCodFreq[38]*sqrtCodFreq[42]*kappa*omega;
	Qij[38*NUM_COD+54] = sqrtCodFreq[38]*sqrtCodFreq[54]*kappa*omega;
	Qij[39*NUM_COD+40] = sqrtCodFreq[39]*sqrtCodFreq[40]*kappa;
	Qij[39*NUM_COD+43] = sqrtCodFreq[39]*sqrtCodFreq[43]*kappa*omega;
	Qij[39*NUM_COD+55] = sqrtCodFreq[39]*sqrtCodFreq[55]*kappa*omega;
	Qij[40*NUM_COD+44] = sqrtCodFreq[40]*sqrtCodFreq[44]*kappa*omega;
	Qij[40*NUM_COD+56] = sqrtCodFreq[40]*sqrtCodFreq[56]*kappa*omega;
	Qij[41*NUM_COD+42] = sqrtCodFreq[41]*sqrtCodFreq[42]*kappa;
	Qij[41*NUM_COD+43] = sqrtCodFreq[41]*sqrtCodFreq[43]*omega;
	Qij[41*NUM_COD+44] = sqrtCodFreq[41]*sqrtCodFreq[44]*omega;
	Qij[41*NUM_COD+57] = sqrtCodFreq[41]*sqrtCodFreq[57]*kappa*omega;
	Qij[42*NUM_COD+43] = sqrtCodFreq[42]*sqrtCodFreq[43]*omega;
	Qij[42*NUM_COD+44] = sqrtCodFreq[42]*sqrtCodFreq[44]*omega;
	Qij[42*NUM_COD+58] = sqrtCodFreq[42]*sqrtCodFreq[58]*kappa*omega;
	Qij[43*NUM_COD+44] = sqrtCodFreq[43]*sqrtCodFreq[44]*kappa;
	Qij[43*NUM_COD+59] = sqrtCodFreq[43]*sqrtCodFreq[59]*kappa*omega;
	Qij[44*NUM_COD+60] = sqrtCodFreq[44]*sqrtCodFreq[60]*kappa*omega;
	Qij[45*NUM_COD+46] = sqrtCodFreq[45]*sqrtCodFreq[46]*kappa;
	Qij[45*NUM_COD+47] = sqrtCodFreq[45]*sqrtCodFreq[47];
	Qij[45*NUM_COD+48] = sqrtCodFreq[45]*sqrtCodFreq[48];
	Qij[45*NUM_COD+49] = sqrtCodFreq[45]*sqrtCodFreq[49]*kappa*omega;
	Qij[45*NUM_COD+53] = sqrtCodFreq[45]*sqrtCodFreq[53]*omega;
	Qij[45*NUM_COD+57] = sqrtCodFreq[45]*sqrtCodFreq[57]*omega;
	Qij[46*NUM_COD+47] = sqrtCodFreq[46]*sqrtCodFreq[47];
	Qij[46*NUM_COD+48] = sqrtCodFreq[46]*sqrtCodFreq[48];
	Qij[46*NUM_COD+50] = sqrtCodFreq[46]*sqrtCodFreq[50]*kappa*omega;
	Qij[46*NUM_COD+54] = sqrtCodFreq[46]*sqrtCodFreq[54]*omega;
	Qij[46*NUM_COD+58] = sqrtCodFreq[46]*sqrtCodFreq[58]*omega;
	Qij[47*NUM_COD+48] = sqrtCodFreq[47]*sqrtCodFreq[48]*kappa;
	Qij[47*NUM_COD+51] = sqrtCodFreq[47]*sqrtCodFreq[51]*kappa*omega;
	Qij[47*NUM_COD+55] = sqrtCodFreq[47]*sqrtCodFreq[55]*omega;
	Qij[47*NUM_COD+59] = sqrtCodFreq[47]*sqrtCodFreq[59]*omega;
	Qij[48*NUM_COD+52] = sqrtCodFreq[48]*sqrtCodFreq[52]*kappa*omega;
	Qij[48*NUM_COD+56] = sqrtCodFreq[48]*sqrtCodFreq[56]*omega;
	Qij[48*NUM_COD+60] = sqrtCodFreq[48]*sqrtCodFreq[60]*omega;
	Qij[49*NUM_COD+50] = sqrtCodFreq[49]*sqrtCodFreq[50]*kappa;
	Qij[49*NUM_COD+51] = sqrtCodFreq[49]*sqrtCodFreq[51];
	Qij[49*NUM_COD+52] = sqrtCodFreq[49]*sqrtCodFreq[52];
	Qij[49*NUM_COD+53] = sqrtCodFreq[49]*sqrtCodFreq[53]*omega;
	Qij[49*NUM_COD+57] = sqrtCodFreq[49]*sqrtCodFreq[57]*omega;
	Qij[50*NUM_COD+51] = sqrtCodFreq[50]*sqrtCodFreq[51];
	Qij[50*NUM_COD+52] = sqrtCodFreq[50]*sqrtCodFreq[52];
	Qij[50*NUM_COD+54] = sqrtCodFreq[50]*sqrtCodFreq[54]*omega;
	Qij[50*NUM_COD+58] = sqrtCodFreq[50]*sqrtCodFreq[58]*omega;
	Qij[51*NUM_COD+52] = sqrtCodFreq[51]*sqrtCodFreq[52]*kappa;
	Qij[51*NUM_COD+55] = sqrtCodFreq[51]*sqrtCodFreq[55]*omega;
	Qij[51*NUM_COD+59] = sqrtCodFreq[51]*sqrtCodFreq[59]*omega;
	Qij[52*NUM_COD+56] = sqrtCodFreq[52]*sqrtCodFreq[56]*omega;
	Qij[52*NUM_COD+60] = sqrtCodFreq[52]*sqrtCodFreq[60]*omega;
	Qij[53*NUM_COD+54] = sqrtCodFreq[53]*sqrtCodFreq[54]*kappa;
	Qij[53*NUM_COD+55] = sqrtCodFreq[53]*sqrtCodFreq[55]*omega;
	Qij[53*NUM_COD+56] = sqrtCodFreq[53]*sqrtCodFreq[56]*omega;
	Qij[53*NUM_COD+57] = sqrtCodFreq[53]*sqrtCodFreq[57]*kappa*omega;
	Qij[54*NUM_COD+55] = sqrtCodFreq[54]*sqrtCodFreq[55]*omega;
	Qij[54*NUM_COD+56] = sqrtCodFreq[54]*sqrtCodFreq[56]*omega;
	Qij[54*NUM_COD+58] = sqrtCodFreq[54]*sqrtCodFreq[58]*kappa*omega;
	Qij[55*NUM_COD+56] = sqrtCodFreq[55]*sqrtCodFreq[56]*kappa;
	Qij[55*NUM_COD+59] = sqrtCodFreq[55]*sqrtCodFreq[59]*kappa*omega;
	Qij[56*NUM_COD+60] = sqrtCodFreq[56]*sqrtCodFreq[60]*kappa*omega;
	Qij[57*NUM_COD+58] = sqrtCodFreq[57]*sqrtCodFreq[58]*kappa;
	Qij[57*NUM_COD+59] = sqrtCodFreq[57]*sqrtCodFreq[59];
	Qij[57*NUM_COD+60] = sqrtCodFreq[57]*sqrtCodFreq[60];
	Qij[58*NUM_COD+59] = sqrtCodFreq[58]*sqrtCodFreq[59];
	Qij[58*NUM_COD+60] = sqrtCodFreq[58]*sqrtCodFreq[60];
	Qij[59*NUM_COD+60] = sqrtCodFreq[59]*sqrtCodFreq[60]*kappa;

	/* Fill in the lower triangle */
	for(i=0;i<NUM_COD;i++) for(j=i+1;j<NUM_COD;j++) Qij[j*NUM_COD+i] = Qij[i*NUM_COD+j];
}

/**
 * Ensures that frequencies are not smaller than MINFREQ and
 * that two frequencies differ by at least 2*MINFDIFF.
 * This avoids potential problems later when eigenvalues
 * are computed.
 */
void CheckCodFrequencies()
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
	for (i = 0; i < NUM_COD; i++) {
		if (codFreq[i] < MINFREQ) codFreq[i] = MINFREQ;
		if (codFreq[i] > maxfreq) {
			maxfreq = codFreq[i];
			maxi = i;
		}
		sum += codFreq[i];
	}
	
	diff = 1.0 - sum;
	codFreq[maxi] += diff;

	for (i = 0; i < NUM_COD - 1; i++) {
		for (j = i+1; j < NUM_COD; j++) {
			if (codFreq[i] == codFreq[j]) {
				codFreq[i] += MINFDIFF;
				codFreq[j] -= MINFDIFF;
			}
		}
	}
	
	for(i=0; i<NUM_COD; i++) sqrtCodFreq[i] = sqrt(codFreq[i]);
}

void SetCodMatrix(double *matrix, double len)
{	
	int i,j,k;
	double expt[NUM_COD];
	double *P;
	
/* P(t)ij = SUM Cijk * exp{Root*t}
*/
	P=matrix;
	if (len<1e-6) { 
		for (i=0; i<NUM_COD; i++) {
			for (j=0; j<NUM_COD; j++) {
				if (i==j)
					*P=1.0;
				else 	
					*P=0.0;
				P++;
			}
		}
		return; 
	}
	
	for (k=1; k<NUM_COD; k++) {
		expt[k]=exp(len*Root[k]);
	}
	for (i=0; i<NUM_COD; i++) {
		double psum = 0.0;
		double *maxP = P;
		for (j=0; j<NUM_COD; j++) {
			(*P)=Cijk[i*SQNUM_COD+j*NUM_COD];
			for (k=1; k<NUM_COD; k++) {
				(*P)+=Cijk[i*SQNUM_COD+j*NUM_COD+k]*expt[k];
			}
			if((*P)>(*maxP)) maxP = P;
			psum += *P;
			P++;
		}
		/* Adjust for rounding error */
		(*maxP) += (1.0-psum);
	}
	
	P=matrix;
	for (i=0; i<NUM_COD; i++) {
		P++;
		for (j=1; j<NUM_COD; j++) {
			*P += *(P - 1);
			P++;
		}
	}
}

void SetCodVector(double *vector, short state, double len)
{	
	int i,j,k;
	double expt[NUM_COD];
	double *P;

	P=vector;
	if (len<1e-6) { 
		for (i=0; i<NUM_COD; i++) {
			if (i==state)
				*P=1.0;
			else 	
				*P=0.0;
			P++;
		}
		return; 
	}
	for (k=1; k<NUM_COD; k++) {
		expt[k]=exp(len*Root[k]);
	}
	
	for (j=0; j<NUM_COD; j++) {
		(*P)=Cijk[state*SQNUM_COD+j*NUM_COD];
		for (k=1; k<NUM_COD; k++) {
			(*P)+=Cijk[state*SQNUM_COD+j*NUM_COD+k]*expt[k];
		}
		P++;
	}
	
	P=vector;
	P++;
	for (j=1; j<NUM_COD; j++) {
		*P += *(P - 1);
		P++;
	}
}


