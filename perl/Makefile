PERL_INCLUDES=/usr/lib64/perl5/CORE

all: perl

###########################################################################################

src/perlforcealigner_wrap.cxx: perlforcealigner.i
	mkdir -p src build
	swig -perl -c++ -I../src -outdir src -o src/perlforcealigner_wrap.cxx perlforcealigner.i
	cp src/perlforcealigner.pm build/.

build/perlforcealigner.so: src/perlforcealigner_wrap.cxx ../src/libgiza.a
	g++  -Wno-deprecated -shared -Xlinker -zmuldefs -fPIC \
	src/perlforcealigner_wrap.cxx \
	../src/libgiza.a \
	-I${PERL_INCLUDES} \
	-I../src \
	-o build/perlforcealigner.so
	
perl: build/perlforcealigner.so

###########################################################################################

clean:
	rm -rf src build