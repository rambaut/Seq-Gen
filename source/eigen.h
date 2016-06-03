/*  Header file for eigen.c                                */

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

#ifndef _EIGEN_H_
#define _EIGEN_H_

int abyx (double a, double x[], int n);
int xtoy (double x[], double y[], int n);
int matinv( double x[], int n, int m, double space[]);
int eigen(int job, double A[], int n, double rr[], double ri[],
          double vr[], double vi[], double w[]);

#endif /* _EIGEN_H_ */
