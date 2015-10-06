/*  Header file for gamma.c                                */

/*  
   Sequence Generator - seq-gen, version 1.3.2
   Andrew Rambaut & Nick Grassly
   Department of Zoology, University of Oxford			
	
   The code in this file is taken from Ziheng Yang's PAML package.
   http://abacus.gene.ucl.ac.uk/

   Any feedback is very welcome.
   http://evolve.zoo.ox.ac.uk/software/Seq-Gen/
   email: andrew.rambaut@zoo.ox.ac.uk
*/

#ifndef _GAMMA_H_
#define _GAMMA_H_

double rndgamma (double s);
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);

#endif /* _GAMMA_H_ */

