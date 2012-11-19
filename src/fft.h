/*---------------------------------------------------------------------------*\

  FILE........: fft.h
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

#include "comp.h"
#ifdef KISS_FFT
#include "kiss_fft.h"
#endif
#ifdef LIBAVCODEC_FFT
#include <libavcodec/avfft.h>
#endif

#ifndef FFT_H_
#define FFT_H_

#ifdef KISS_FFT
typedef kiss_fft_cfg fft_cfg;
#elif defined(LIBAVCODEC_FFT)
typedef struct {
    FFTContext *context;
    int size;
} *fft_cfg;
#else
#error FFT engine was not defined
#endif

fft_cfg fft_new(int n, int inverse);

void fft_delete(const fft_cfg cfg);

void fft_do(const fft_cfg cfg, const COMP *in, COMP *out);

#endif /* FFT_H_ */
