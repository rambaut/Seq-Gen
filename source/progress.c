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
#include "progress.h"

static int barLength, dotGap, dots, bar;

/*********************/
void InitProgressBar(int total)
{
	barLength=total;
	dotGap=1;
	if (barLength>MAX_BAR_LENGTH) {
		dotGap=barLength/MAX_BAR_LENGTH;
		if (barLength%dotGap)
			barLength=(barLength/dotGap)+1;
		else
			barLength/=dotGap;
	}
}


/*********************/
void DrawProgressBar()
{
	int i;

	if (barLength<2)
		return;
		
	fprintf(stderr, "0%%|");
	for (i=0; i<(barLength); i++)
		fputc('_', stderr);
	fprintf(stderr, "|100%%\n  [");
	fflush(stderr);
	dots=0;
	bar=0;
}


/*********************/
void ProgressBar()
{
	if (barLength<2)
		return;

	if (bar%dotGap==0) {
		fputc('.', stderr);
		fflush(stderr);
		dots++;
		if (dots==barLength) {
			fputc(']', stderr);
			fputc('\n', stderr);
		}
		fflush(stderr);
	}
		
	bar++;
}


