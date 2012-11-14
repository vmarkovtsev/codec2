/*--------------------------------------------------------------------------*\

	FILE........: vqtrain_clipped.c
	AUTHOR......: David Rowe
	DATE CREATED: 22 March 2011

	This program trains vector quantisers using K dimensional Lloyd-Max
	method, modified to centre clip errors.

\*--------------------------------------------------------------------------*/

/*
  Copyright (C) 2009 David Rowe

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

/*-----------------------------------------------------------------------*\

				INCLUDES

\*-----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*-----------------------------------------------------------------------*\

				DEFINES

\*-----------------------------------------------------------------------*/

#define	DELTAQ 	0.001		/* quiting distortion			*/
#define	MAX_STR	80		/* maximum string length		*/

/*-----------------------------------------------------------------------*\

			FUNCTION PROTOTYPES

\*-----------------------------------------------------------------------*/

void zero(float v[], int k);
void acc(float v1[], float v2[], int k);
void norm(float v[], int k, long n);
long quantise(float cb[], float vec[], int k, int m, float clip, float *se);

/*-----------------------------------------------------------------------*\

				MAIN

\*-----------------------------------------------------------------------*/

int main(int argc, char *argv[]) {
    long   k,m;		/* dimension and codebook size			*/
    float  *vec;	/* current vector 				*/
    float  *cb;		/* vector codebook				*/
    float  *cent;	/* centroids for each codebook entry		*/
    long   *n;		/* number of vectors in this interval		*/
    long   J;		/* number of vectors in training set		*/
    long   i,j;
    long   ind;	     	/* index of current vector			*/
    float  se;		/* squared error for this iteration		*/
    float  Dn,Dn_1;	/* current and previous iterations distortion	*/
    float  delta;	/* improvement in distortion 			*/
    FILE   *ftrain;	/* file containing training set			*/
    FILE   *fvq;	/* file containing vector quantiser		*/
    float  clip;

    /* Interpret command line arguments */

    if (argc != 6)	{
	printf("usage: vqtrain TrainFile K M VQFile clip\n");
	exit(0);
    }

    /* Open training file */

    ftrain = fopen(argv[1],"rb");
    if (ftrain == NULL) {
	printf("Error opening training database file: %s\n",argv[1]);
	exit(1);
    }

    /* determine k and m, and allocate arrays */

    k = atol(argv[2]);
    m = atol(argv[3]);
    clip = atof(argv[5]);
    printf("dimension K=%ld  number of entries M=%ld clip = %f\n", k,m,clip);
    vec = (float*)malloc(sizeof(float)*k);
    cb = (float*)malloc(sizeof(float)*k*m);
    cent = (float*)malloc(sizeof(float)*k*m);
    n = (long*)malloc(sizeof(long)*m);
    if (cb == NULL || cb == NULL || cent == NULL || vec == NULL) {
	printf("Error in malloc.\n");
	exit(1);
    }

    /* determine size of training set */

    J = 0;
    while(fread(vec, sizeof(float), k, ftrain) == k)
    J++;
    printf("J=%ld entries in training set\n", J);

    /* set up initial codebook state from samples of training set */

    rewind(ftrain);
    fread(cb, sizeof(float), k*m, ftrain);

    /* main loop */

    Dn = 1E32;
    j = 1;
    do {
	Dn_1 = Dn;

	/* zero centroids */

	for(i=0; i<m; i++) {
	    zero(&cent[i*k], k);
	    n[i] = 0;
	}

	/* quantise training set */

	se = 0.0;
	rewind(ftrain);
	for(i=0; i<J; i++) {
	    fread(vec, sizeof(float), k, ftrain);
	    ind = quantise(cb, vec, k, m, clip, &se);
	    n[ind]++;
	    acc(&cent[ind*k], vec, k);
	}
	Dn = se/J;
	delta = (Dn_1-Dn)/Dn;

	printf("\r  Iteration %ld, Dn = %f, Delta = %e\n", j, Dn, delta);
	j++;

	/* determine new codebook from centriods */

	if (delta > DELTAQ)
	    for(i=0; i<m; i++) {
		if (n[i] != 0) {
		    norm(&cent[i*k], k, n[i]);
		    memcpy(&cb[i*k], &cent[i*k], k*sizeof(float));
		}
	    }

    } while (delta > DELTAQ);

    /* save codebook to disk */

    fvq = fopen(argv[4],"wt");
    if (fvq == NULL) {
	printf("Error opening VQ file: %s\n",argv[4]);
	exit(1);
    }

    fprintf(fvq,"%ld %ld\n",k,m);
    for(j=0; j<m; j++) {
	for(i=0; i<k; i++)
	    fprintf(fvq,"%f  ",cb[j*k+i]);
	fprintf(fvq,"\n");
    }
    fclose(fvq);

    return 0;
}

/*-----------------------------------------------------------------------*\

				FUNCTIONS

\*-----------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\

	FUNCTION....: zero()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Zeros a vector of length k.

\*---------------------------------------------------------------------------*/

void zero(float v[], int k)
/*  float  v[];		ptr to start of vector		*/
/*  int    k;		lngth of vector			*/
{
    int	i;

    for(i=0; i<k; i++)
	v[i] = 0.0;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: acc()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Adds k dimensional vectors v1 to v2 and stores the result back in v1.

\*---------------------------------------------------------------------------*/

void acc(float v1[], float v2[], int k)
/*  float  v1[];	ptr to start of vector to accumulate	*/
/*  float  v2[];	ptr to start of vector to add		*/
/*  int	   k;		dimension of vectors			*/
{
    int	   i;

    for(i=0; i<k; i++)
	v1[i] += v2[i];
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: norm()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Divides each element in k dimensional vector v by n.

\*---------------------------------------------------------------------------*/

void norm(float v[], int k, long n)
/*  float  v[];		ptr to start of vector		*/
/*  int	   k;		dimension of vectors		*/
/*  int	   n;		normalising factor		*/
{
    int	   i;

    for(i=0; i<k; i++)
	v[i] /= n;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: quantise()

	AUTHOR......: David Rowe
	DATE CREATED: 23/2/95

	Quantises vec by choosing the nearest vector in codebook cb, and
	returns the vector index.  The squared error of the quantised vector
	is added to se.

        This quantiser uses a centre clipping method.  The quant error
	of each vector element is inspected.  If it is beneath a
	certain threshold it is ignored in the error measurement.

\*---------------------------------------------------------------------------*/

long quantise(float cb[], float vec[], int k, int m, float clip, float *se)
/* float   cb[][K];	current VQ codebook		*/
/* float   vec[];	vector to quantise		*/
/* int	   k;		dimension of vectors		*/
/* int     m;		size of codebook		*/
/* float   clip;        centre clip level               */
/* float   *se;		accumulated squared error 	*/
{
    float   element_e;
    float   e;		   /* error for current vector	*/
    long    besti;         /* best index so far		*/
    float   beste;	   /* best error so far		*/
    long    j;
    int     i;

    besti = 0;
    beste = 1E32;
    for(j=0; j<m; j++) {
	e = 0.0;
	for(i=0; i<k; i++) {
	    element_e = cb[j*k+i]-vec[i];
	    if (fabs(element_e) > clip) {
		e +=  element_e*element_e;
	    }
	}
	if (e < beste) {
	    beste = e;
	    besti = j;
	}
    }

    *se += beste;

    return(besti);
}

