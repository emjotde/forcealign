#include <iostream>
#include <fstream>

#include "SymForceAligner.h"

int main(int argc, char* argv[]) {
    std::string srcLang = argv[2];
    std::string trgLang = argv[3];
    std::string mosesModelPath = argv[1];
    size_t CPUS = 20;
    
    SymForceAligner fa(srcLang, trgLang, mosesModelPath, CPUS);
    fa.setMode(SymForceAligner::GrowDiagFinal);
    
    ifstream srcIn;
    ifstream trgIn;
    
    srcIn.open(argv[4]);
    trgIn.open(argv[5]);
    
    std::string srcLine;
    std::string trgLine;
    std::cerr << "Adding sentences" << std::endl;
    while(std::getline(srcIn, srcLine) && std::getline(trgIn, trgLine))
        fa.addSentence(srcLine, trgLine);
	std::cerr << "Aligning corpus" << std::endl;
    fa.alignCorpus();
    
    ofstream out;
    out.open(argv[6]);
    out << fa.getAlignmentsStr();
    
	return 0;
}