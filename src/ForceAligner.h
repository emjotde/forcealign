#ifndef FORCEALIGNER_H__
#define FORCEALIGNER_H__

#include "getSentence.h"
#include "TTables.h"
#include "model1.h"
#include "model2.h"
#include "model3.h"
#include "defs.h"
#include "vocab.h"
#include "Perplexity.h"
#include "Dictionary.h"
#include "Parameter.h"
#include "myassert.h"
#include "D4Tables.h"
#include "D5Tables.h"

#include <set>
#include <vector>
#include <sstream>
#include <map>

typedef std::pair<size_t, size_t> AlignmentPoint;
typedef std::set<AlignmentPoint>  Alignment;
typedef std::vector<Alignment>    Alignments;

typedef std::pair<WordIndex, WordIndex> Cooccurrence;
typedef std::set<Cooccurrence>          Cooccurrences;


class SentenceHandlerMem : public sentenceHandler {
  private:
    std::vector<std::string> m_e;
    std::vector<std::string> m_f;
    
    std::vector<sentPair> m_buffer;
    
    vcbList* m_elist;
    vcbList* m_flist;
    
    size_t m_rewinds;
    
    Cooccurrences cooccurrences;
    
    int getNextSentence(sentPair& sent, vcbList* elist, vcbList* flist);
    WordIndex getId(std::string word, vcbList* vocab);
    //bool readNextSentence(sentPair& sent);
        
  public:
    SentenceHandlerMem(vcbList* elist=0, vcbList* flist=0);
    virtual Cooccurrences& getCooccurrences() { return cooccurrences; }
    void rewind();
    void addSentence(std::string e, std::string f);
    

};

class ForceAligner {
  private:
    std::string m_eVcbFile;
    std::string m_fVcbFile;
    std::string m_coocFile;
    std::string m_prefix;
    
    SentenceHandlerMem m_corpus;    
    
    typedef std::pair<WordIndex, WordIndex> WordPair;
    
    vcbList* m_eVcbList;
    vcbList* m_fVcbList;
    Dictionary* m_dictionary;
    std::map<WordPair, PROB>* m_tTable;
    amodel<PROB>* m_aTable;
    amodel<COUNT>* m_aCountTable;
    nmodel<PROB>* m_nTable;
    amodel<PROB>* m_dTable;
    
    WordClasses* m_eWordClasses;
    WordClasses* m_fWordClasses;
    
    d4model* m_d4m;
    
    double m_p0;
    double m_p1;
    
    void init();
    bool readProbTable(std::string &filename);
    void filterTTable(tmodel<COUNT, PROB> &tempTTable, Cooccurrences &cooc);
    
  public:
    static Mutex s_mutex;
    
    static Alignments s_alignments;
    Alignments m_alignments;
    
    //ForceAligner(std::string eVcbFile, std::string fVcbFile, std::string prefix);
    ForceAligner(std::string src, std::string trg, std::string path, size_t CPUS = 8);
    
    ~ForceAligner();
    
    void addSentence(std::string e, std::string f);
    void alignCorpus();
    
    Alignment& alignSentence(std::string e, std::string f);
    std::string alignSentenceStr(std::string e, std::string f);
    
    void clearAll();
    
    Alignments& getAlignments();
    std::string getAlignmentsStr();
};

#endif
