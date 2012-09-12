#ifndef SYMFORCEALIGNER_H__
#define SYMFORCEALIGNER_H__

#include <set>
#include <vector>

#include "ForceAligner.h"

#define MAX_WORD 1000  //maximum lengthsource/target strings 
#define MAX_M 200     //maximum length of source strings
#define MAX_N 200     //maximum length of target strings

#define VERBOSE 0

class SymForceAligner {
  private:
    bool m_diagonal;
    bool m_final;    
    bool m_bothuncovered;
    
    Alignments m_alignments;
    
    ForceAligner m_srcTrgAligner;
    ForceAligner m_trgSrcAligner;
    
    Alignment cUnion(int *a, int m, int* b, int n);
    Alignment cIntersection(int *a, int m, int* b, int n);
    Alignment cGrow(int *a, int m, int *b, int n);
    Alignment symmetrize(Alignment& a_, Alignment& b_);
    
  public:
    enum Mode { Src2Trg, Trg2Src, Intersection, Union, Grow, GrowDiag, GrowDiagFinal, GrowDiagFinalAnd };
    
  private:
    Mode m_mode;
    
  public:
    SymForceAligner(std::string src, std::string trg, std::string path, size_t CPUS = 8);
    void setMode(Mode mode);
    
    void addSentence(std::string e, std::string f);
    void alignCorpus();
    
    Alignment& alignSentence(std::string e, std::string f);
    std::string alignSentenceStr(std::string e, std::string f);
    
    void clearAll();
    
    Alignments& getAlignments();
    std::string getAlignmentsStr();
};

#endif
