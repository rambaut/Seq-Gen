/*  Header file for gamma.c                                */

/*  
   Sequence Generator - seq-gen, version 1.3.5
   Copyright (c)1996-2025, Andrew Rambaut
   Institute of Evolutionary Biology, University of Edinburgh			

   Any feedback is very welcome.
   http://tree.bio.ed.ac.uk/software/seqgen/
   email: a.rambaut@ed.ac.uk
*/

#ifndef _GAMMA_H_
#define _GAMMA_H_

double rndgamma (double s);
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);

#endif /* _GAMMA_H_ */

