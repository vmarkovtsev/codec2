/*---------------------------------------------------------------------------*\
                                                 
  FILE........: fft.c                                                  
  AUTHOR......: Bruce Robertson                                      
  DATE CREATED: 20/11/2010                            
                                                         
  Bridging function to the kiss_fft package.      
                                                               
\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2010 Bruce Robertson

  All rights reserved.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License version 2.1, as
  published by the Free Software Foundation.  This program is
  distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

/*
 Added support for libavcodec FFT - Markovtsev Vadim <v.markovtsev@samsung.com>
 Added some minor NEON optimizations - Markovtsev Vadim <v.markovtsev@samsung.com>
*/


#include <assert.h>
#ifdef NEON
#include <arm_neon.h>
#endif
#ifdef KISS_FFT
#include <kiss_fft.h>

/*---------------------------------------------------------------------------*\
                                                                            
                                GLOBALS                                       
                                                                             
\*---------------------------------------------------------------------------*/

kiss_fft_cpx *fin;
kiss_fft_cpx *fout;
kiss_fft_cfg cfg_forward;
kiss_fft_cfg cfg_reverse;

/*---------------------------------------------------------------------------*\
                                                                             
  initialize_fft(int n)                                                                  
                                                                             
  Initialisation function for kiss_fft. This assumes that all calls to fft() 
  use the same datatypes and are one arrays of the same size.

\*---------------------------------------------------------------------------*/

void
initialize_fft (int n)
{
  fin = KISS_FFT_MALLOC (n * sizeof (kiss_fft_cpx));
  assert(fin != NULL);
  fout = KISS_FFT_MALLOC (n * sizeof (kiss_fft_cpx));
  assert(fout != NULL);
  cfg_forward = kiss_fft_alloc (n, 0, NULL, NULL);
  assert(cfg_forward != NULL);
  cfg_reverse = kiss_fft_alloc (n, 1, NULL, NULL);
  assert(cfg_reverse != NULL);
}

/*---------------------------------------------------------------------------*\
                                                                             
  fft(float x[], int n, int isign)                                                
  Function that calls kiss_fft with the signature of four1 from NRC.

\*---------------------------------------------------------------------------*/


void
fft (float x[], int n, int isign)
{
  if (cfg_forward == NULL)
    {
      initialize_fft (n);
    }
  int c;
#ifndef NEON
  for (c = 0; c < n * 2; c += 2)
    {
      fin[c / 2].r = x[c];
      fin[c / 2].i = -x[c + 1];
    }
#else
  const float32_t conjmem[4] = { 1.0f, -1.0f, 1.0f, -1.0f };
  const float32x4_t conjvec = vld1q_f32(conjmem);
  for (c = 0; c < n; c += 2)
		{
			float32x4_t cmplxpair = vld1q_f32(&x[c * 2]);
			cmplxpair = vmulq_f32(cmplxpair, conjvec);
			vst1q_f32((float32_t *)&fin[c], cmplxpair);
		}
#endif
  kiss_fft_cfg cfg;
  if (isign == -1)
    {
      cfg = cfg_reverse;
    }
  else
    {
      cfg = cfg_forward;
    }
  kiss_fft (cfg, fin, fout);
#ifndef NEON
  for (c = 0; c < n * 2; c += 2)
    {
      x[c] = fout[(c) / 2].r;
      x[c + 1] = -fout[(c) / 2].i;
    }
#else
  for (c = 0; c < n; c += 2)
		{
			float32x4_t cmplxpair = vld1q_f32((float32_t *)&fout[c]);
			cmplxpair = vmulq_f32(cmplxpair, conjvec);
			vst1q_f32(&x[c * 2], cmplxpair);
		}
#endif
}

#else  // #ifdef KISS_FFT
#ifdef LIBAVCODEC_FFT
#include <malloc.h>
#include <libavcodec/avfft.h>

FFTComplex *fftbuffer = NULL;
FFTContext *fftcontext[2] = { NULL, NULL };

int log2int(int n)
{
	int res = 0;
	while (n >>= 1)
	  {
	    res++;
	  }
	return res;
}

FFTContext *
initialize_fft (int n, int isign)
{
	static int length = 0;
	int inverse = (isign == -1);
	if (length < n || fftbuffer == NULL)
	  {
		  if (fftbuffer == NULL)
		    {
	         fftbuffer = malloc(sizeof(FFTComplex) * n);
		    }
		  else
		    {
			    fftbuffer = realloc(fftbuffer, sizeof(FFTComplex) * n);
		    }
		  length = n;
		  if (fftcontext[inverse] != NULL)
				{
					av_fft_end(fftcontext[inverse]);
					fftcontext[inverse] = NULL;
				}
	  }
	if (fftcontext[inverse] == NULL)
		{
			fftcontext[inverse] = av_fft_init(log2int(n), inverse);
			assert(fftcontext[inverse] != NULL && "av_fft_init() failed");
		}
	return fftcontext[inverse];
}

void
fft (float x[], int n, int isign)
{
	FFTContext *cntxt = initialize_fft (n, isign);
	int c;
#ifndef NEON
	for (c = 0; c < n * 2; c += 2)
		{
			fftbuffer[c / 2].re = x[c];
			fftbuffer[c / 2].im = -x[c + 1];
		}
#else
	const float32_t conjmem[4] = { 1.0f, -1.0f, 1.0f, -1.0f };
	const float32x4_t conjvec = vld1q_f32(conjmem);
	for (c = 0; c < n; c += 2)
		{
			float32x4_t cmplxpair = vld1q_f32(&x[c * 2]);
			cmplxpair = vmulq_f32(cmplxpair, conjvec);
			vst1q_f32((float32_t *)&fftbuffer[c], cmplxpair);
		}
#endif
	av_fft_permute(cntxt, fftbuffer);
	av_fft_calc(cntxt, fftbuffer);
#ifndef NEON
	for (c = 0; c < n * 2; c += 2)
		{
			x[c] = fftbuffer[c / 2].re;
			x[c + 1] = -fftbuffer[c / 2].im;
		}
#else
	for (c = 0; c < n; c += 2)
		{
			float32x4_t cmplxpair = vld1q_f32((float32_t *)&fftbuffer[c]);
			cmplxpair = vmulq_f32(cmplxpair, conjvec);
			vst1q_f32(&x[c * 2], cmplxpair);
		}
#endif
}

#else

#error No FFT implementation was chosen!

#endif  // #ifdef LIBAVCODEC_FFT
#endif  // #ifdef KISS_FFT
