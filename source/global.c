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
#include <string.h>
#include <ctype.h>

#include "global.h"

int verbose=0, verboseMemory=0, quiet=0, userSeed;
long totalMem=0;
unsigned long randomSeed;

/* functions */

/*************************************/
void *AllocMem(long n, char *name, char *func, int showInfo)
{
	void *P;
	
	if ( (P=malloc(n))==NULL ) {
		fprintf(stderr, "Out of memory allocating '%s': %s()\n", name, func);
		exit(0);
	}
	
	totalMem+=n;

	if (showInfo && verboseMemory)
		fprintf(stderr, "%s in %s() - %ld bytes\n", name, func, n);
	
	return P;
}


/*************************************/
void *CAllocMem(long n, char *name, char *func, int showInfo)
{
	void *P;
	
	if ( (P=calloc(n, 1))==NULL ) {
		fprintf(stderr, "Out of memory allocating '%s': %s()\n", name, func);
		exit(0);
	}
	
	totalMem+=n;
	if (showInfo && verboseMemory)
		fprintf(stderr, "%s in %s() - %ld bytes\n", name, func, n);
	
	return P;
}


/*************************************/
int GetDoubleParams(int argc, char **argv, int *argn, char *pos, int numParams, double *params)
{
	int i;
	char *st, buf[256];
	
	i=0;
	strcpy(buf, pos);
	st=strtok(buf, "\t,/");
	do {
		if (st==NULL) {
			if ((*argn)+1<argc) {
				(*argn)++;
				strcpy(buf, argv[*argn]);
			} else
				return -1;
			st=strtok(buf, "\t,/");
			if (st==NULL)
				return -1;
		}
		
		if (sscanf(st, "%lf", params+i)!=1)
			return -1;
		i++;
		if (i<numParams)
			st=strtok(NULL, "\t,/");
	} while (i<numParams);

	return 0;
}


/*************************************/
int GetIntParams(int argc, char **argv, int *argn, char *pos, int numParams, int *params)
{
	int i;
	char *st, buf[256];
	
	i=0;
	strcpy(buf, pos);
	st=strtok(buf, "\t,/");
	do {
		if (st==NULL) {
			if ((*argn)+1<argc) {
				(*argn)++;
				strcpy(buf, argv[*argn]);
			} else
				return -1;
			st=strtok(buf, "\t,/");
			if (st==NULL)
				return -1;
		}
		
		if (sscanf(st, "%d", params+i)!=1)
			return -1;
		i++;
		if (i<numParams)
			st=strtok(NULL, "\t,/");
	} while (i<numParams);

	return 0;
}


/*************************************/
int GetUnsignedLongParams(int argc, char **argv, int *argn, char *pos, int numParams, unsigned long *params)
{
	int i;
	char *st, buf[256];
	
	i=0;
	strcpy(buf, pos);
	st=strtok(buf, "\t,/");
	do {
		if (st==NULL) {
			if ((*argn)+1<argc) {
				(*argn)++;
				strcpy(buf, argv[*argn]);
			} else
				return -1;
			st=strtok(buf, "\t,/");
			if (st==NULL)
				return -1;
		}
		
		if (sscanf(st, "%lu", params+i)!=1)
			return -1;
		i++;
		if (i<numParams)
			st=strtok(NULL, "\t,/");
	} while (i<numParams);

	return 0;
}


/*************************************/
int GetStrParam(int argc, char **argv, int *argn, char *pos, char *param, int len)
{
	int i;
	char *st, *P, buf[256];
	
	i=0;
	strcpy(buf, pos);
	st=strtok(buf, "\t,/");
	if (st==NULL) {
		if ((*argn)+1<argc) {
			(*argn)++;
			strcpy(buf, argv[*argn]);
		} else
			return -1;
		st=strtok(buf, "\t,/");
		if (st==NULL)
			return -1;
	}
	strncpy(param, st, len);
	param[len]='\0';
	P=param;
	while (*P) {
		*P=toupper(*P);
		P++;
	}
	
	return 0;
}
