/*---------------------------------------------------------------------------*\

  FILE........: fft.c
  AUTHOR......: Markovtsev Vadim
  DATE CREATED: 19/11/12

  FFT engine abstraction layer.

\*---------------------------------------------------------------------------*/

/*
  Copyright (C) 2009 David Rowe

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

#include "fft.h"
#include <assert.h>
#ifdef LIBAVCODEC_FFT
#include <malloc.h>
#include <string.h>
#endif

#ifdef LIBAVCODEC_FFT
static int log2int(int n)
{
    int res = 0;
    while (n >>= 1)
      {
        res++;
      }
    return res;
}
#endif

fft_cfg fft_new(int n, int inverse) {
#ifdef KISS_FFT
    kiss_fft_cfg cfg = kiss_fft_alloc (n, inverse, NULL, NULL);
    return cfg;
#elif defined(LIBAVCODEC_FFT)
    FFTContext *ctxt = av_fft_init(log2int(n), inverse);
    if (ctxt == NULL) return NULL;
    fft_cfg cfg = malloc(sizeof(*cfg));
    cfg->context = ctxt;
    cfg->size = sizeof(COMP) * n;
    return cfg;
#else
#error FFT engine was not defined
#endif
}

void fft_delete(const fft_cfg cfg) {
#ifdef KISS_FFT
    KISS_FFT_FREE(cfg);
#elif defined(LIBAVCODEC_FFT)
    av_fft_end(cfg->context);
    free(cfg);
#else
#error FFT engine was not defined
#endif
}

void fft_do(const fft_cfg cfg, const COMP *in, COMP *out) {
#ifdef KISS_FFT
    fft_do(cfg, in, out);
#elif defined(LIBAVCODEC_FFT)
    memcpy(out, in, cfg->size);
    av_fft_permute(cfg->context, (FFTComplex *)out);
    av_fft_calc(cfg->context, (FFTComplex *)out);
#else
#error FFT engine was not defined
#endif
}
