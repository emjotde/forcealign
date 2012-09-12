JAVA_INCLUDES=/usr/lib/jvm/java-1.6.0-sun-1.6.0.30.x86_64/include

JAVA_SOURCES=example.java src/javaforcealignerJNI.java src/ForceAligner.java src/SymForceAligner.java
JAVA_CLASSES=build/example.class

all: java

###########################################################################################

src/javaforcealigner_wrap.cxx: javaforcealigner.i
	mkdir -p src build
	swig -java -c++ -I../src -outdir src -o src/javaforcealigner_wrap.cxx javaforcealigner.i

build/libjavaforcealigner.so: src/javaforcealigner_wrap.cxx ../src/libgiza.a
	g++ -Wno-deprecated -shared -Xlinker -zmuldefs -fPIC \
	src/javaforcealigner_wrap.cxx \
	../src/libgiza.a \
	-I${JAVA_INCLUDES} -I${JAVA_INCLUDES}/linux \
	-I../src \
	-o build/libjavaforcealigner.so
	
${JAVA_CLASSES}: ${JAVA_SOURCES}
	javac -d build ${JAVA_SOURCES}
	
java: build/libjavaforcealigner.so ${JAVA_CLASSES}

###########################################################################################

clean:
	rm -rf src build