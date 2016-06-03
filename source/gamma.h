/*  Header file for gamma.c                                */

/*  
   Sequence Generator - seq-gen, version 1.3.3
   Andrew Rambaut & Nick Grassly
   Institute of Evolutionary Biology, University of Edinburgh			

   The code in this file is taken from Ziheng Yang's PAML package.
   http://abacus.gene.ucl.ac.uk/

   Any feedback is very welcome.
   http://tree.bio.ed.ac.uk/software/seqgen/
   email: a.rambaut@ed.ac.uk
*/

#ifndef _GAMMA_H_
#define _GAMMA_H_

double rndgamma (double s);
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);

#endif /* _GAMMA_H_ */

