if LIBAVCODEC_FFT
FFTCFLAGS=@LIBAVCODEC_CFLAGS@
FFTLIBS=@LIBAVCODEC_LIBS@
endif
if MATH_NEON
MATHNEONCFLAGS=-I$(top_srcdir)/math-neon/include
MATHNEONLIBS=$(top_builddir)/math-neon/libmath_neon.la
endif

AM_CFLAGS = -I$(top_srcdir)/src -fPIC -g  -DFLOATING_POINT -DVAR_ARRAYS -O2 -Wall $(MATHNEONCFLAGS) $(FFTCFLAGS)
AUTOMAKE_OPTS = gnu
NAME = libcodec2
AM_CPPFLAGS = $(AM_CFLAGS)

noinst_PROGRAMS = genres genlsp extract vqtrain vqtrainjnd tnlp tinterp tquant vq_train_jvm scalarlsptest tfdmdv t48_8 lspsync create_interleaver tlspsens vqtrainph genphdata genampdata polar2rect vqtrainsp tprede pre de 

genres_SOURCES = genres.c ../src/lpc.c
genres_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
genres_LDFLAGS = $(LIBS)

genlsp_SOURCES = genlsp.c ../src/lpc.c ../src/lsp.c
genlsp_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
genlsp_LDFLAGS = $(LIBS)

extract_SOURCES = extract.c
extract_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
extract_LDFLAGS = $(LIBS)

vqtrain_SOURCES = vqtrain.c
vqtrain_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
vqtrain_LDFLAGS = $(LIBS)

vqtrainjnd_SOURCES = vqtrainjnd.c
vqtrainjnd_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
vqtrainjnd_LDFLAGS = $(LIBS)

vqtrainph_SOURCES = vqtrainph.c
vqtrainph_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
vqtrainph_LDFLAGS = $(LIBS)

vqtrainsp_SOURCES = vqtrainsp.c
vqtrainsp_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
vqtrainsp_LDFLAGS = $(LIBS)

genphdata_SOURCES = genphdata.c
genphdata_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
genphdata_LDFLAGS = $(LIBS)

genampdata_SOURCES = genampdata.c
genampdata_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
genampdata_LDFLAGS = $(LIBS)

polar2rect_SOURCES = polar2rect.c
polar2rect_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
polar2rect_LDFLAGS = $(LIBS)

vq_train_jvm_SOURCES = vq_train_jvm.c
vq_train_jvm_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
vq_train_jvm_LDFLAGS = $(LIBS)

CODEBOOKS = ../src/codebook.c ../src/codebookd.c ../src/codebookvq.c ../src/codebookjnd.c ../src/codebookdt.c ../src/codebookjvm.c ../src/codebookvqanssi.c ../src/codebookge.c

tnlp_SOURCES = tnlp.c ../src/sine.c ../src/nlp.c ../src/kiss_fft.c ../src/fft.c ../src/dump.c ../src/comp.c
tnlp_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
tnlp_LDFLAGS = $(LIBS)

tinterp_SOURCES = tinterp.c ../src/sine.c ../src/kiss_fft.c ../src/fft.c ../src/interp.c ../src/lpc.c ../src/lsp.c ../src/quantise.c $(CODEBOOKS) ../src/dump.c ../src/comp.c
tinterp_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
tinterp_LDFLAGS = $(LIBS)

tquant_SOURCES = tquant.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/comp.c ../src/kiss_fft.c ../src/fft.c $(CODEBOOKS)
tquant_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
tquant_LDFLAGS = $(LIBS)

scalarlsptest_SOURCES = scalarlsptest.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/comp.c ../src/kiss_fft.c ../src/fft.c $(CODEBOOKS)
scalarlsptest_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
scalarlsptest_LDFLAGS = $(LIBS)

tfdmdv_SOURCES = tfdmdv.c ../src/fdmdv.c ../src/kiss_fft.c ../src/fft.c ../src/octave.c ../src/comp.c
tfdmdv_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
tfdmdv_LDFLAGS = $(LIBS)

t48_8_SOURCES = t48_8.c ../src/fdmdv.c ../src/kiss_fft.c ../src/fft.c ../src/comp.c
t48_8_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
t48_8_LDFLAGS = $(LIBS)

lspsync_SOURCES = lspsync.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/comp.c ../src/kiss_fft.c ../src/fft.c \
../src/codec2.c ../src/sine.c ../src/nlp.c ../src/postfilter.c ../src/phase.c ../src/interp.c ../src/pack.c $(CODEBOOKS)
lspsync_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
lspsync_LDFLAGS = $(LIBS)

create_interleaver_SOURCES = create_interleaver.c 
create_interleaver_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
create_interleaver_LDFLAGS = $(LIBS)

tlspsens_SOURCES = tlspsens.c ../src/quantise.c ../src/lpc.c ../src/lsp.c ../src/dump.c ../src/comp.c ../src/kiss_fft.c ../src/fft.c ../src/codec2.c ../src/sine.c ../src/nlp.c ../src/pack.c ../src/interp.c ../src/postfilter.c ../src/phase.c $(CODEBOOKS)
tlspsens_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
tlspsens_LDFLAGS = $(LIBS)

tprede_SOURCES = tprede.c ../src/lpc.c 
tprede_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
tprede_LDFLAGS = $(LIBS)

pre_SOURCES = pre.c ../src/lpc.c 
pre_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
pre_LDFLAGS = $(LIBS)

de_SOURCES = de.c ../src/lpc.c 
de_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS)
de_LDFLAGS = $(LIBS)

#tsrc_SOURCES = tsrc.c 
#tsrc_LDADD = $(lib_LTLIBRARIES) $(MATHNEONLIBS) $(FFTLIBS) /usr/lib/libsamplerate.so.0
#tsrc_LDFLAGS = $(LIBS)
