/*
 *  codonmodels.h
 *  
 *
 *  Created by Daniel Wilson on 3/6/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _CODONMODELS_H_
#define _CODONMODELS_H_

#include "evolve.h"

#define NUM_COD 61
#define SQNUM_COD 3721
#define CUNUM_COD 226981

enum { COD_NONE = -1, COD_NY98, numCodModels };

extern char *codons;

extern int codFreqSet;

extern double codFreq[NUM_COD];
extern double codAddFreq[NUM_COD];
extern double codMatrix[MAX_RATE_CATS][SQNUM_COD];
extern double codVector[NUM_COD];

extern double kappa, omega;

void SetCodModel(int theModel);
void SetCodMatrix(double *matrix, double len);
void SetCodVector(double *vector, short state, double len);

#endif //_CODONMODELS_H_
