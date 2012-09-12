## Process this file with automake to produce Makefile.in

## Created by Anjuta

INCLUDES = \
	-DPACKAGE_LOCALE_DIR=\""$(prefix)/$(DATADIRNAME)/locale"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(datadir)"\" 

AM_CFLAGS =\
	 -Wall -Wno-deprecated

bin_PROGRAMS = forcealign

forcealign_SOURCES = \
	forcealign.cpp

forcealign_DEPENDENCIES = \
	libgiza.a

forcealign_CXXFLAGS = \
	-MT -Wno-deprecated -Xlinker -zmuldefs \
	-MD \
	-MP \
	-MF \
	-fPIC

forcealign_CPPFLAGS = \
	-DNDEBUG \
	-DWORDINDEX_WITH_4_BYTE \
	-DBINARY_SEARCH_FOR_TTABLE \
	-DDEBUG

forcealign_LDFLAGS = \
	-L.

forcealign_LDADD = \
	-lgiza \
	-lpthread

lib_LIBRARIES = \
	libgiza.a

libgiza_a_SOURCES = \
	alignment.cpp\
	alignment.h \
	AlignTables.cpp \
	AlignTables.h \
	Array.h \
	Array2.h \
	Array4.h \
	ATables.cpp \
	ATables.h \
	collCounts.cpp \
	collCounts.h \
	common.h \
	D4Tables.h \
	D5Tables.h \
	defs.h \
	Dictionary.cpp \
	Dictionary.h \
	file_spec.h \
	FlexArray.h \
	ForceAligner.h \
	ForceAligner.cpp \
	SymForceAligner.h \
	SymForceAligner.cpp \
	ForwardBackward.cpp \
	ForwardBackward.h \
	getSentence.cpp \
	getSentence.h \
	Globals.h \
	hmm.cpp \
	hmm.h \
	HMMTables.cpp \
	HMMTables.h \
	logprob.cpp \
	logprob.h \
	model1.cpp \
	model1.h \
	model2.cpp \
	model2.h \
	model2to3.cpp \
	model3.cpp \
	model3.h \
	model3_viterbi.cpp \
	model3_viterbi_with_tricks.cpp \
	model345-peg.cpp \
	MoveSwapMatrix.cpp \
	MoveSwapMatrix.h \
	myassert.cpp \
	myassert.h \
	mymath.h \
	mystl.h \
	NTables.cpp \
	NTables.h \
	Parameter.cpp \
	Parameter.h \
	parse.cpp \
	Perplexity.cpp \
	Perplexity.h \
	Pointer.h \
	reports.cpp \
	SetArray.cpp \
	SetArray.h \
	syncObj.h \
	transpair_model1.h \
	transpair_model2.h \
	transpair_model3.cpp \
	transpair_model3.h \
	transpair_model4.cpp \
	transpair_model4.h \
	transpair_model5.cpp \
	transpair_model5.h \
	transpair_modelhmm.h \
	ttableDiff.hpp \
	TTables.cpp \
	TTables.h \
	types.h \
	utility.cpp \
	utility.h \
	Vector.h \
	vocab.cpp \
	vocab.h \
	WordClasses.h

libgiza_a_CXXFLAGS = \
	-MD -Wno-deprecated \
	-MP \
	-MF \
	-MT \
	-fPIC

libgiza_a_CPPFLAGS = \
	-DNDEBUG \
	-DWORDINDEX_WITH_4_BYTE \
	-DBINARY_SEARCH_FOR_TTABLE \
	-DDEBUG