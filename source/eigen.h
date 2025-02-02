/*  Header file for eigen.c                                */

/*  
   Sequence Generator - seq-gen, version 1.3.5
   Copyright (c)1996-2025, Andrew Rambaut
   Institute of Evolutionary Biology, University of Edinburgh			

   The code in this file is taken from Ziheng Yang's PAML package, distributed under the GNU GPL v3
   http://abacus.gene.ucl.ac.uk/
   https://github.com/abacus-gene/paml

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
