/*  Header file for aamodels.c                             */

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

#ifndef _AA_MODELS_H_
#define _AA_MODELS_H_

#include "evolve.h"

#define NUM_AA 20
#define SQNUM_AA 400
#define CUNUM_AA 8000
#define NUM_AA_REL_RATES 190

enum { AA_NONE = -1, AA_JTT, AA_WAG, AA_DAYHOFF78, AA_BLOSUM62, AA_MTREV24, AA_CPREV45, AA_MTART, AA_GENERAL, numAAModels };

extern char *aminoAcids;

enum { ala, arg, asn, asp, cys, gln, glu, gly, his, ileu, leu, lys, met, phe, pro, ser, thr, trp, tyr, val};

extern int aaFreqSet;

extern double aaFreq[NUM_AA];
extern double aaAddFreq[NUM_AA];
extern double aaMatrix[MAX_RATE_CATS][SQNUM_AA];
extern double aaVector[NUM_AA];

extern double aaRelativeRate[NUM_AA_REL_RATES];

void SetAAModel(int theAAModel);
void SetAAMatrix(double *matrix, double len);
void SetAAVector(double *vector, short state, double len);

#endif /* _AA_MODELS_H_ */
