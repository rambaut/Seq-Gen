/*  Header file to define tree and node structures         */

/*  
   Sequence Generator - seq-gen, version 1.3.2
   Copyright (c)1996-2004, Andrew Rambaut & Nick Grassly
   Department of Zoology, University of Oxford			
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
   http://evolve.zoo.ox.ac.uk/software/Seq-Gen/
   email: andrew.rambaut@zoo.ox.ac.uk
*/
 
#ifndef _TREE_H_
#define _TREE_H_

#define MAX_NAME_LEN 256

typedef struct TNode TNode;
struct TNode {
	TNode *branch0, *branch1, *branch2, *next;
	double length0, length1, length2, param;
	int tipNo;
	
	char *sequence;
};


typedef struct TTree TTree;
struct TTree {
	int rooted, lengths;
	TNode *root, *nodeList;
	int numTips, numNodes;
	double totalLength;
	char **names;
	TNode **tips;
	int capacity;
};

#endif /* _TREE_H_ */

