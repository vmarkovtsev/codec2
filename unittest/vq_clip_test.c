/*--------------------------------------------------------------------------*\

	FILE........: vq_clip_test.c
	AUTHOR......: David Rowe
	DATE CREATED: 22 March 2011

	This program tests a vector quantiser against a test database 
	using the centre cliped quantisation method.

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

#define THRESH	40.0	/* threshold energy/sample for frame inclusion 	*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defines.h"
#include "comp.h"
#include "codec2_internal.h"	
#include "codec2.h"	
#include "interp.h"	
#include "lpc.h"	

long quantise_clip(float cb[], float vec[], int k, int m, float clip, float *se,
	      int *less_than_clip);

float sample_log_amp_quad_nl(
			     float w[],     /* frequency points            */
			     float A[],     /* for these amplitude samples */
			     int np,        /* number of frequency points  */
			     float w_sample /* frequency of new samples    */
			     );

int switch_present(sw,argc,argv)
  char sw[];     /* switch in string form */
  int argc;      /* number of command line arguments */
  char *argv[];  /* array of command line arguments in string form */
{
  int i;       /* loop variable */

  for(i=1; i<argc; i++)
    if (!strcmp(sw,argv[i]))
      return(i);

  return 0;
}

struct codebook {
  unsigned int	k;
  unsigned int	log2m;
  unsigned int	m;
  float * cb;
};

float
get_float(FILE * in, const char * name, char * * cursor, char * buffer,
 int size)
{
  for ( ; ; ) {
    char *	s = *cursor;
    char	c;

    while ( (c = *s) != '\0' && !isdigit(c) && c != '-' && c != '.' )
      s++;
     
    /* Comments start with "#" and continue to the end of the line. */
    if ( c != '\0' && c != '#' ) {
      char *	end = 0;
      float	f = 0;

      f = strtod(s, &end);

      if ( end != s )
        *cursor = end;
        return f;
    }

    if ( fgets(buffer, size, in) == NULL ) {
      fprintf(stderr, "%s: Format error.\n", name);
      exit(1);
    }
    *cursor = buffer;
  }
}

static struct codebook *
load(FILE * file, const char * name)
{
  char			line[1024];
  char *		cursor = line;
  struct codebook *	b = malloc(sizeof(struct codebook));
  int			i;
  int			size;

  *cursor = '\0';

  b->k = (int)get_float(file, name, &cursor, line, sizeof(line));
  b->m = (int)get_float(file, name ,&cursor, line, sizeof(line));
  size = b->k * b->m;

  b->cb = (float *)malloc(size * sizeof(float));

  for ( i = 0; i < size; i++ )
    b->cb[i] = get_float(file, name, &cursor, line, sizeof(line));

  return b;
}

/*---------------------------------------------------------------------------*\
                                                       
  FUNCTION....: synthesise_one_frame_orig_phase()	     
  AUTHOR......: David Rowe			      
  DATE CREATED: 23/8/2010 

  Synthesise 80 speech samples (10ms) from model parameters.

\*---------------------------------------------------------------------------*/

void synthesise_one_frame_orig_phase(CODEC2 *c2, short speech[], MODEL *model)
{
    int     i;

    postfilter(model, &c2->bg_est);
    synthesise(c2->Sn_, model, c2->Pn, 1);

    for(i=0; i<N; i++) {
	if (c2->Sn_[i] > 32767.0)
	    speech[i] = 32767;
	else if (c2->Sn_[i] < -32767.0)
	    speech[i] = -32767;
	else
	    speech[i] = c2->Sn_[i];
    }
}

#define PHASE_POINTS 50

int main(int argc, char *argv[]) {
    FILE   *fraw, *fvq, *fout;
    short   buf[N];
    void   *c2;
    CODEC2 *c3;
    MODEL   model;
    float w[MAX_AMP], A[MAX_AMP];
    float wres[MAX_AMP], Ares[MAX_AMP], AresdB[MAX_AMP], AresdB_q[MAX_AMP], AresdB_prev[MAX_AMP];
    float deltat[MAX_AMP];
    float sam, E;
    long  f, af;
    int   i, j, err_less_than_clip_low, err_less_than_clip_high, k, m, besti_low, besti_high;
    float av_err_low, av_err_high;
    float clip;
    struct  codebook *cblow = malloc(sizeof(struct codebook));
    struct  codebook *cbhigh = malloc(sizeof(struct codebook));
    float sum_av_err_low, sum_av_err_high;
    int  sum_err_less_than_clip_low,sum_err_less_than_clip_high;
    float signal, noise, snr, new_A, sum_snr;
    float  wph[PHASE_POINTS], Pph[PHASE_POINTS], R[LPC_ORD+1], aks[LPC_ORD+1];
    int phase_model;

    if (argc < 5) {
	printf("usage: vq_clip_test RawFileIn LowVQFile HighVQFile clip [RawFileOut] --phase\n");
	exit(0);
    }

    f = af = 0;

    c2 = codec2_create();
    c3 = (CODEC2*)c2;

    for(i=0; i<RES_POINTS; i++)
	AresdB_prev[i] = 0.0; 

    /* Open files */

    fraw = fopen(argv[1],"rb");
    if (fraw == NULL) {
	printf("Error opening input .raw file: %s",argv[1]);
	exit(1);
    }

    fvq = fopen(argv[2],"rt");
    if (fvq == NULL) {
	printf("Error opening input VQ file: %s",argv[2]);
	exit(1);
    }
    cblow = load(fvq, argv[2]);
    fclose(fvq);

    fvq = fopen(argv[3],"rt");
    if (fvq == NULL) {
	printf("Error opening input VQ file: %s",argv[3]);
	exit(1);
    }
    cbhigh = load(fvq, argv[3]);
    fclose(fvq);

    clip = atof(argv[4]);

    if (argc >= 6) {
	fout = fopen(argv[5],"wb");
	if (fout == NULL) {
	    printf("Error opening output .raw file: %s",argv[4]);
	    exit(1);
	}
    }
    else
	fout = NULL;

    if (switch_present("--phase",argc, argv)) {
	phase_model = 1;
    }
    else
	phase_model = 0;

#ifdef DUMP
    dump_on("test");
#endif

    sum_av_err_low = sum_av_err_high = 0;
    sum_err_less_than_clip_low = sum_err_less_than_clip_high = 0;
    sum_snr = 0;

    while(fread(buf, sizeof(short), N, fraw) == N) {
	analyse_one_frame(c3, &model, buf);

	E = 0.0;
	for(i=0; i<N; i++) {
	    sam = (float)buf[i];
	    E += sam*sam;
	}
	E = 10.0*log10(E/N);

	f++;
	printf("Frame: %ld  \n", f);

	resample_amp_fixed(&model, w, A, wres, Ares,
			   AresdB_prev, AresdB, deltat);
#ifdef DUMP
	dump_resample(wres,Ares,RES_POINTS);
#endif
	//#define VQ
#ifdef VQ
	besti_low = quantise_clip(cblow->cb, deltat, cblow->k, cblow->m, clip, 
				  &av_err_low, 
				  &err_less_than_clip_low);
	printf("LOW av_err: %3.2f err_less_than_clip: %d  besti: %d\n\n", 
	       av_err_low, err_less_than_clip_low, besti_low );

	besti_high = quantise_clip(cbhigh->cb, &deltat[cblow->k], cbhigh->k, cbhigh->m, clip, 
				  &av_err_high, 
				  &err_less_than_clip_high);
	printf("HIGH av_err: %3.2f err_less_than_clip: %d  besti: %d\n\n", 
	       av_err_high, err_less_than_clip_high, besti_high );

	/* recover Ares vector */

	for(i=0; i<cblow->k; i++) {
	    AresdB_q[i] = AresdB_prev[i] + cblow->cb[besti_low*cblow->k+i];
	    Ares[i] = pow(10.0, AresdB_q[i]/20.0);
	    //printf("%d %f 2%f\n", i, AresdB[i], AresdB_q[i]);
	}
	for(i=0,j=cblow->k; i<cbhigh->k; i++,j++) {
	    AresdB_q[j] = AresdB_prev[j] + cbhigh->cb[besti_high*cbhigh->k+i];
	    Ares[j] = pow(10.0, AresdB_q[j]/20.0);
	    // printf("%d %f %f\n", j, AresdB[j], AresdB_q[j]);
	}
#endif

	//#define ART_NOISE
#ifdef ART_NOISE
	// add some artificial noise to test effect on phase model

	for(i=0; i<RES_POINTS; i++) {
	    noise = 3.0*(1.0 - 2.0*(float)rand()/RAND_MAX);
	    AresdB[i] = 20.0*log10(Ares[i]) + noise;
	    Ares[i] = pow(10.0, AresdB[i]/20.0);
	}
#endif

	/* recover {A} vector */

	signal = noise = 0.0;
    
	for(i=1; i<=model.L; i++) {
	    //printf("model.Wo*i %f\n", model.Wo*i);
	    new_A = pow(10.0,sample_log_amp_quad_nl(wres, Ares, RES_POINTS, model.Wo*i));
	    signal += pow(model.A[i], 2.0);
	    noise  += pow(model.A[i] - new_A, 2.0);
	    model.A[i] = new_A;
	}

	snr = 10.0*log10(signal/noise);
	printf("snr = %3.2f\n", snr);

	if (phase_model) {
	    /* generate a phase spectrum using freq domain LPC */

	    for(i=0; i<PHASE_POINTS; i++) {
		wph[i] = i*PI/PHASE_POINTS;
		new_A = pow(10.0,sample_log_amp_quad_nl(wres, Ares, RES_POINTS, wph[i]));
		Pph[i] = new_A*new_A;
	    }

	    //for(i=0; i<RES_POINTS; i++)
	    //	printf("%d wres: %f Ares: %f Pw %f\n", i, wres[i], Ares[i], Pph[i]);

	    {
		COMP    Sw[FFT_ENC];
		dft_speech(Sw, c3->Sn, c3->w);

#ifdef DUMP
		dump_Sn(c3->Sn); dump_Sw(Sw); dump_model(&model);
#endif
	    }

	    //	    #define TIME_DOM
#ifdef TIME_DOM
	    {
		float Wn[M];
		for(i=0; i<M; i++)
		    Wn[i] = c3->Sn[i]*c3->w[i];
		autocorrelate(Wn,R,M,LPC_ORD);
	    }
	    levinson_durbin(R, aks, LPC_ORD);
#ifdef DUMP
	    dump_ak(aks, LPC_ORD);
#endif
#else
	    autocorrelate_freq(Pph, wph, R, PHASE_POINTS, LPC_ORD);
	    levinson_durbin(R, aks, LPC_ORD);
#ifdef DUMP
	    dump_ak(aks, LPC_ORD);
#endif
#endif

	    for(i=0; i<=LPC_ORD; i++)
		printf("%d R: %f aks: %f\n", i, R[i], aks[i]);
	    //if (f==30)
	    //	exit(0);
	    synthesise_one_frame(c3, buf, &model, aks);
	}
	else {
	    synthesise_one_frame_orig_phase(c3, buf, &model);
	}

	if (fout != NULL) fwrite(buf, sizeof(short), N, fout);
	
	/* If energy high enough, include this frame in stats */

	if (E > THRESH) {
	    af++;
	    sum_av_err_low += av_err_low;
	    if (err_less_than_clip_low)
		sum_err_less_than_clip_low++;
	    sum_av_err_high += av_err_high;
	    if (err_less_than_clip_high)
		sum_err_less_than_clip_high++;
	    sum_snr += snr;
	}
	for(i=0; i<RES_POINTS; i++)
	    AresdB_prev[i] = AresdB_q[i];	    
    }

    printf("Active frames: %ld av_error_low: %3.2f  %% less than clip low %3.2f %%"
	   "av_error_high: %3.2f  %% less than clip high %3.2f %%SNR av %3.2f\n", 
	   af, 
	   sum_av_err_low/af, 
	   100.0*(float)sum_err_less_than_clip_low/af,
	   sum_av_err_high/af, 
	   100.0*(float)sum_err_less_than_clip_high/af,
	   sum_snr/af);

    codec2_destroy(c2);

    fclose(fraw);
    if (fout != NULL) fclose(fout);

    return 0;
}

/*---------------------------------------------------------------------------*\

	FUNCTION....: quantise_clip()

	AUTHOR......: David Rowe
	DATE CREATED: 23 March 2011

	Quantises vec by choosing the nearest vector in codebook cb, and
	returns the vector index.  

        This quantiser uses a centre clipping method.  The quant error
	of each vector element is inspected.  If it is beneath a
	certain threshold it is ignored in the error measurement.

\*---------------------------------------------------------------------------*/

long quantise_clip(float cb[], float vec[], int k, int m, float clip, float *av_err,
	      int *less_than_clip)
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
		e += fabs(element_e);
	    }
	}
	if (e < beste) {
	    beste = e;
	    besti = j;
	}
    }

    *av_err = 0;
    *less_than_clip = 1;
    for(i=0; i<k; i++) {
	element_e = cb[besti*k+i]-vec[i];
	(*av_err) += fabs(element_e);
	if (fabs(element_e) > clip) {
	    *less_than_clip = 0;
	}
    }
    for(i=0; i<k; i++)
	printf("% 14f", vec[i]);
    printf("\n");
    for(i=0; i<k; i++)
	printf("% 14f", cb[besti*k+i]);
    printf("\n");

    (*av_err) /= k;

    return(besti);
}
