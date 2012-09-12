/* File: perlforcealigner.i */
%include "std_string.i"

%module perlforcealigner
%{
#include "ForceAligner.h"
#include "SymForceAligner.h"
%}

class ForceAligner {
  public:
    ForceAligner(std::string src, std::string trg, std::string path);
    
    void addSentence(std::string e, std::string f);
    void alignCorpus();
    std::string alignSentenceStr(std::string e, std::string f);
    void clearAll();
    std::string getAlignmentsStr();
};

class SymForceAligner {
  public:
    enum Mode { Src2Trg, Trg2Src, Intersection, Union, Grow, GrowDiag, GrowDiagFinal, GrowDiagFinalAnd };
    SymForceAligner(std::string src, std::string trg, std::string path);
    void setMode(Mode mode);
    
    void addSentence(std::string e, std::string f);
    void alignCorpus();
    std::string alignSentenceStr(std::string e, std::string f);
    void clearAll();
    std::string getAlignmentsStr();
};
