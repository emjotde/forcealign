#include "ForceAligner.h"

extern short NCPUS;

Alignments ForceAligner::s_alignments;
Mutex ForceAligner::s_mutex;

void ForceAligner::init() {
    Vector<WordEntry> *evlist, *fvlist;
    evlist = new Vector<WordEntry>();
    fvlist = new Vector<WordEntry>();
    
    m_eVcbList = new vcbList(*evlist);
    m_fVcbList = new vcbList(*fvlist);

    std::cerr << "reading vocabulary files \n";
    
    m_eVcbList->setName(m_eVcbFile.c_str());
    m_fVcbList->setName(m_fVcbFile.c_str());
    m_eVcbList->readVocabList();
    m_fVcbList->readVocabList();       
    
    std::cerr << "Source vocabulary list has " << m_eVcbList->uniqTokens()
                    << " unique tokens \n";

    std::cerr << "Target vocabulary list has " << m_fVcbList->uniqTokens()
                    << " unique tokens \n";

    std::string prev_t    = m_prefix + ".t3.final";
    std::string prev_p0   = m_prefix + ".p0_3.final";
    std::string prev_a    = m_prefix + ".a3.final";
    std::string prev_d    = m_prefix + ".d3.final";
    std::string prev_n    = m_prefix + ".n3.final";

    m_tTable = new std::map<WordPair, PROB>();
    std::cerr << "We are going to load previous model " << prev_t << endl;
    if(!readProbTable(prev_t)){
        cerr << "Failed reading " << prev_t << endl; 
        exit(1);
    }   

    m_aTable = new amodel<PROB>(false);
    m_aCountTable = new amodel<COUNT>(false);
    
    cerr << "We are going to load previous model from " << prev_a << endl;
    if(!m_aTable->readTable(prev_a.c_str())){
            cerr << "Failed reading " << prev_a << endl;
            exit(1);
    }

    m_nTable = new nmodel<PROB>(m_eVcbList->size()+1, MAX_FERTILITY);
    cerr << "We are going to load previous N model from " << prev_n << endl;
    if(!m_nTable->readNTable(prev_n.c_str())){
            cerr << "Failed reading " << prev_n << endl;
            exit(1);
    }
    
    m_dTable = new amodel<PROB>(true);
    cerr << "We are going to load previous D model from " << prev_d << endl;
    if(!m_dTable->readTable(prev_d.c_str())){
            cerr << "Failed reading " << prev_d << endl;
            exit(1);
    }
    
    cerr << "We are going to load previous P0 Value model from " << prev_p0 << endl;
    ifstream ifs(prev_p0.c_str());
    ifs >> m_p0;
    m_p1 = 1 - m_p0;
    
    std::string prev_d4   = m_prefix + ".d4.final";
    std::string prev_d4_2 = m_prefix + ".D4.final";
    
    m_eWordClasses = new WordClasses();
    m_fWordClasses = new WordClasses();
    
    m_d4m = new d4model(MAX_SENTENCE_LENGTH, *m_eWordClasses, *m_fWordClasses);
    
    std::cerr << "We are going to read d4 table from " << prev_d4 << "," << prev_d4_2  << std::endl;
    m_d4m->readProbTable(prev_d4.c_str(), prev_d4_2.c_str());
    
    m_d4m->makeWordClasses(*m_eVcbList, *m_fVcbList, m_eVcbFile + ".classes",
                    m_fVcbFile + ".classes", *m_eVcbList, *m_fVcbList);
    
    m_corpus = SentenceHandlerMem(m_eVcbList, m_fVcbList);
}

bool ForceAligner::readProbTable(std::string &filename) {
    ifstream inf(filename.c_str());
    cerr << "Reading t prob. table from " << filename << "\n";
    if (!inf) {
        std::cerr << "\nERROR: Cannot open " << filename << "\n";
        return false;
    }
    WordIndex src_id, trg_id;
    PROB prob;
    int nEntry=0;
    while (inf >> src_id >> trg_id >> prob) {
            (*m_tTable)[WordPair(src_id, trg_id)] = prob;
            nEntry++;
    }
    cerr << "Read " << nEntry << " entries in prob. table.\n";
    return true;
}

//ForceAligner::ForceAligner(std::string eVcbFile, std::string fVcbFile,
//             std::string coocFile, std::string prefix)
//: m_eVcbFile(eVcbFile), m_fVcbFile(fVcbFile), m_coocFile(coocFile),
//  m_prefix(prefix), m_eVcbList(0), m_fVcbList(0), m_dictionary(0),
//  m_tTable(0), m_aTable(0), m_aCountTable(0), m_nTable(0), m_dTable(0),
//  m_p0(1), m_p1(0), m_eWordClasses(0), m_fWordClasses(0), m_d4m(0)
//{
//    init();
//}

ForceAligner::ForceAligner(std::string src, std::string trg, std::string path, size_t CPUS)
: m_eVcbList(0), m_fVcbList(0), m_corpus(m_eVcbList, m_fVcbList), m_dictionary(0),
  m_tTable(0), m_aTable(0), m_aCountTable(0), m_nTable(0), m_dTable(0),
  m_p0(1), m_p1(0), m_eWordClasses(0), m_fWordClasses(0), m_d4m(0)
{
    NCPUS = CPUS;
    
    std::string suf = trg + "-" + src;
    
    m_eVcbFile = path + "/corpus/" + src + ".vcb";
    m_fVcbFile = path + "/corpus/" + trg + ".vcb";
    m_coocFile = path + "/giza." + suf + "/" + suf + ".cooc";
    m_prefix   = path + "/giza." + suf + "/" + suf;
    
    init();
}

ForceAligner::~ForceAligner() {
    delete m_eVcbList;
    delete m_fVcbList;
    delete m_dictionary;
    delete m_tTable;
    delete m_aTable;
    delete m_aCountTable;
    delete m_nTable;
    delete m_dTable;
}

void ForceAligner::filterTTable(tmodel<COUNT, PROB> &tempTTable, Cooccurrences &cooc) {
    typedef tmodel<COUNT, PROB>::CPPair CPPair;
        
    tempTTable.recordDiff = false;
    
    std::pair<unsigned int,CPPair> cp;
    std::vector<std::pair<unsigned int,CPPair> > cps;
    int count = 0, count2 = 0;
    int e, f, olde = -1, oldf = -1;
    for(Cooccurrences::iterator it = cooc.begin(); it != cooc.end(); it++) {
        WordIndex e = it->first;
        WordIndex f = it->second;
        
        cp.first = f;
        assert(e >= olde);
        assert(e > olde || f > oldf);
        if(e != olde && olde >= 0 ){
            int oldsize = tempTTable.lexmat.size();
            tempTTable.lexmat.resize(olde + 1);
            for(unsigned int i = oldsize; i < tempTTable.lexmat.size(); ++i)
                tempTTable.lexmat[i] = 0;
            tempTTable.lexmat[olde] = new std::vector<std::pair<unsigned int,CPPair> > (cps);
            
            cps.clear();
            count2 += tempTTable.lexmat[olde]->capacity();
        }
        cps.push_back(cp);
        olde = e;
        oldf = f;
        count++;
    }
    tempTTable.lexmat.resize(olde + 1);
    tempTTable.lexmat[olde] = new std::vector<std::pair<unsigned int,CPPair> >(cps);
    count2 += tempTTable.lexmat[olde]->capacity();
    
    tempTTable.mutex.resize(tempTTable.lexmat.size());
    
    for(Cooccurrences::iterator it = cooc.begin(); it != cooc.end(); it++) {
        WordIndex e = it->first;
        WordIndex f = it->second;
        
        tempTTable.insert(e, f, 0, (*m_tTable)[WordPair(e, f)]);
    }
};  

void ForceAligner::addSentence(std::string e, std::string f) {
    m_corpus.addSentence(e, f);
}

void ForceAligner::clearAll() {
    s_alignments.clear();
    m_alignments.clear();
}

void ForceAligner::alignCorpus() {
    s_alignments.clear();
    
    std::cerr << "Reading in corpus and updating vocabulary list" << std::endl;
    sentenceHandler* testCorpus = 0;
    
    std::cerr << "Filtering and updating tables" << std::endl;
    tmodel<COUNT, PROB> tempTTable;
    filterTTable(tempTTable, m_corpus.getCooccurrences());
    // update n_table !
    
    std::cerr << "Preparing models" << std::endl;
    Perplexity trainPerp, testPerp, trainViterbiPerp, testViterbiPerp;
    model1 m1("dummy", *m_eVcbList, *m_fVcbList, tempTTable,
                    trainPerp, m_corpus, &testPerp, testCorpus, trainViterbiPerp,
                    &testViterbiPerp);
    
    model2 m2(m1, *m_aTable, *m_aCountTable);
    model3 m3(m2, *m_dTable, *m_nTable);
    m3.p0 = m_p0;
    m3.p1 = m_p1;
    
    d5model d5m(*m_d4m);
    
    std::string align = "none";
    std::string modelname = "Model4";
    char from = '4';
    char to   = '4';
    
    m_corpus.rewind();
    
    //m3.viterbi_thread(0, align, true, *m_d4m, d5m, true, from, to, modelname);    
    
    std::cerr << "Preparing and running threads" << std::endl;
    vector<m3_em_loop_t> th;
    th.resize(NCPUS);
    
    for(size_t k = 0; k < NCPUS ; k++) {
        th[k].m = &m3;
    	th[k].d4 = m_d4m;
    	th[k].d5 = &d5m;
        th[k].done = 0;
        th[k].valid = 0;
        th[k].it = 0;
        th[k].final = true;
        th[k].alignfile = align;
        th[k].dump_files = true;
        th[k].fromModel = from;
        th[k].toModel = to;
        th[k].modelName = modelname;
        
        th[k].valid = pthread_create(&(th[k].thread), NULL, m3_exe_emloop, &(th[k]));
        
        if(th[k].valid)
            cerr << "Error starting thread " << k << endl;
    }
    
    for (size_t k=0; k<NCPUS; k++){
        pthread_join((th[k].thread), NULL);
        cerr << "Thread " << k << "done" << endl;
    }
    
    m_alignments = s_alignments;
    s_alignments.clear();
}

Alignments& ForceAligner::getAlignments() {
    return m_alignments;
}

std::string ForceAligner::getAlignmentsStr() {
    std::stringstream ss;
    for(Alignments::iterator it1 = m_alignments.begin(); it1 != m_alignments.end(); it1++) {
        Alignment a = *it1;
        for(Alignment::iterator it2 = a.begin(); it2 != a.end(); it2++) {
            if(it2 != a.begin())
                ss << " ";
            ss << it2->first << "-" << it2->second;
        }
        ss << std::endl;
    }
    return ss.str();
}

Alignment& ForceAligner::alignSentence(std::string e, std::string f) {
    clearAll();
    addSentence(e, f);
    alignCorpus();
    return getAlignments()[0];
}

std::string ForceAligner::alignSentenceStr(std::string e, std::string f) {
    Alignment a = alignSentence(e, f);
    
    std::stringstream ss;
    for(Alignment::iterator it = a.begin(); it != a.end(); it++) {
        if(it != a.begin())
            ss << " ";
        ss << it->first << "-" << it->second;
    }
    return ss.str();
}

/*********************************************************************************************/

SentenceHandlerMem::SentenceHandlerMem(vcbList* elist, vcbList* flist)
: m_elist(elist), m_flist(flist), m_rewinds(0)
{
    pthread_mutex_init(&readsent_mutex,NULL);
    pthread_mutex_init(&setprob_mutex,NULL);
    
    position = 0;
    readflag = false ;
    allInMemory = false;
    currentSentence = 0;

    totalPairs1 = 0;
    totalPairs2 = 0;
    pair_no = 0 ;
    noSentInBuffer = 0;
    Buffer.clear();
}

void SentenceHandlerMem::addSentence(std::string e, std::string f) {
    sentPair sent;
    sent.noOcc = 1;
    sent.realCount = sent.noOcc;
    
    std::string word;  
    
    istrstream src(e.c_str());
    sent.eSent.push_back(0);
    while(src >> word) {
        WordIndex w = getId(word, m_elist);
        sent.eSent.push_back(w);
    }
    
    istrstream trg(f.c_str());
    sent.fSent.push_back(0);
    while(trg >> word) {
        WordIndex w = getId(word, m_flist);
        sent.fSent.push_back(w);
    }
    
    if ((sent.fSent.size()-1) > (MAX_FERTILITY-1) * (sent.eSent.size()-1)){
        std::cerr << "WARNING: The following sentence pair has source/target sentence length ration more than\n"<<
            "the maximum allowed limit for a source word fertility\n"<<
            " source length = " << sent.eSent.size()-1 << " target length = " << sent.fSent.size()-1 <<
            " ratio " << double(sent.fSent.size()-1)/  (sent.eSent.size()-1) << " ferility limit : " <<
            MAX_FERTILITY-1 << '\n';
        std::cerr << "Shortening sentence \n";
        std::cerr << sent;
        sent.eSent.resize(min(sent.eSent.size(),sent.fSent.size()));
        sent.fSent.resize(min(sent.eSent.size(),sent.fSent.size()));
    }
    
    if (m_elist && m_flist){
        if ((*m_elist).size() > 0)
            for (WordIndex i= 0 ; i < sent.eSent.size() ; i++){
                if (sent.eSent[i] >= (*m_elist).uniqTokens()){
                    std::cerr << "ERROR: source word " << sent.eSent[i] << " is not in the vocabulary list \n";
                    exit(-1);
                }
                (*m_elist).incFreq(sent.eSent[i], sent.realCount);
            }
        if ((*m_flist).size() > 0)
            for (WordIndex j= 1 ; j < sent.fSent.size() ; j++){
                if (sent.fSent[j] >= (*m_flist).uniqTokens()){
                    std::cerr << "ERROR: target word " << sent.fSent[j] << " is not in the vocabulary list \n";
                    exit(-1);
                }
                (*m_flist).incFreq(sent.fSent[j], sent.realCount);
            }
    }
    
    for(std::vector<WordIndex>::iterator it1 = sent.eSent.begin(); it1 != sent.eSent.end(); it1++)
        for(std::vector<WordIndex>::iterator it2 = sent.fSent.begin(); it2 != sent.fSent.end(); it2++)
            if(*it2 > 0)
                cooccurrences.insert(Cooccurrence(*it1, *it2));

    sent.sentenceNo = ++pair_no;
    m_buffer.push_back(sent);
    
    totalPairs1++;
    totalPairs2 += sent.realCount;
    
    if(m_buffer.size() % 10000 == 0)
        std::cerr << "[" << m_buffer.size() << "]" << std::endl;
}

int SentenceHandlerMem::getNextSentence(sentPair& sent, vcbList* elist, vcbList* flist)
{
    pthread_mutex_lock(&readsent_mutex);
    pair_no++;
    if(pair_no <= m_buffer.size()) {
        sent = m_buffer[pair_no-1];
        position++;
        pthread_mutex_unlock(&readsent_mutex);
        return position;
    }
    else {
        sent.eSent.clear();
        sent.fSent.clear();
        sent.sentenceNo = 0 ;
        sent.noOcc = 0 ;
        sent.realCount=0;
        pthread_mutex_unlock(&readsent_mutex);
        return 0;
    }
    
    //pthread_mutex_lock(&readsent_mutex);
    // 
    //sentPair s;
    //if(readNextSentence(s)){
    //    if ((s.fSent.size()-1) > (MAX_FERTILITY-1) * (s.eSent.size()-1)){
    //        std::cerr << "WARNING: The following sentence pair has source/target sentence length ration more than\n"<<
    //            "the maximum allowed limit for a source word fertility\n"<<
    //            " source length = " << s.eSent.size()-1 << " target length = " << s.fSent.size()-1 <<
    //            " ratio " << double(s.fSent.size()-1)/  (s.eSent.size()-1) << " ferility limit : " <<
    //            MAX_FERTILITY-1 << '\n';
    //        std::cerr << "Shortening sentence \n";
    //        std::cerr << s;
    //        s.eSent.resize(min(s.eSent.size(),s.fSent.size()));
    //        s.fSent.resize(min(s.eSent.size(),s.fSent.size()));
    //    }
    //    
    //    if(m_rewinds == 0) {
    //        if (elist && flist){
    //            if ((*elist).size() > 0)
    //                for (WordIndex i= 0 ; i < s.eSent.size() ; i++){
    //                    if (s.eSent[i] >= (*elist).uniqTokens()){
    //                        std::cerr << "ERROR: source word " << s.eSent[i] << " is not in the vocabulary list \n";
    //                        exit(-1);
    //                    }
    //                    (*elist).incFreq(s.eSent[i], s.realCount);
    //                }
    //            if ((*flist).size() > 0)
    //                for (WordIndex j= 1 ; j < s.fSent.size() ; j++){
    //                    if (s.fSent[j] >= (*flist).uniqTokens()){
    //                        std::cerr << "ERROR: target word " << s.fSent[j] << " is not in the vocabulary list \n";
    //                        exit(-1);
    //                    }
    //                    (*flist).incFreq(s.fSent[j], s.realCount);
    //                }
    //        }
    //    }
    //    position++;
    //    
    //    if(m_buffer.size() <= pair_no)
    //        m_buffer.resize(pair_no + 1);    
    //    m_buffer[pair_no] = s;
    //    sent = m_buffer[pair_no];
    //    
    //    pthread_mutex_unlock(&readsent_mutex);
    //    return position ;
    //}
    //pthread_mutex_unlock(&readsent_mutex);
    //return 0;
}

WordIndex SentenceHandlerMem::getId(std::string word, vcbList* vocab) {
    if(!vocab->has_word(word)) {
        Vector<WordEntry> &list = vocab->getVocabList();
        WordIndex id = list.size();
        list.push_back(WordEntry(word, 0));
        return id;
    }
    else
        return (*vocab)(word);
}

//bool SentenceHandlerMem::readNextSentence(sentPair& sent) {
//    if(pair_no < m_e.size() && pair_no < m_f.size()) { 
//        sent.clear();
//        sent.noOcc = 1;
//        sent.realCount = sent.noOcc;
//        
//        std::string word;  
//        
//        istrstream src(m_e[pair_no].c_str());
//        sent.eSent.push_back(0);
//        while(src >> word) {
//            WordIndex w = getId(word, m_elist);
//            sent.eSent.push_back(w);
//        }
//        
//        istrstream trg(m_f[pair_no].c_str());
//        sent.fSent.push_back(0);
//        while(trg >> word) {
//            WordIndex w = getId(word, m_flist);
//            sent.fSent.push_back(w);
//        }
//        
//        if(m_rewinds == 0) {
//            for(std::vector<WordIndex>::iterator it1 = sent.eSent.begin(); it1 != sent.eSent.end(); it1++)
//                for(std::vector<WordIndex>::iterator it2 = sent.fSent.begin(); it2 != sent.fSent.end(); it2++)
//                    if(*it2 > 0)
//                        cooccurrences.insert(Cooccurrence(*it1, *it2));
//        }
//        
//        sent.sentenceNo = ++pair_no;
//        return true;
//    }
//    else {
//        sent.eSent.clear();
//        sent.fSent.clear();
//        sent.sentenceNo = 0 ;
//        sent.noOcc = 0 ;
//        sent.realCount=0;
//        return false;
//    }
//}

void SentenceHandlerMem::rewind() {
    m_rewinds++;
    position = 0;
    currentSentence = 0;
    readflag = false ;
    pair_no = 0 ;
    noSentInBuffer = 0 ;
    Buffer.clear();
}
    
/*********************************************************************************************/

void printAlignToFile(const Vector<WordIndex>& es,
                      const Vector<WordIndex>& fs,
                      const Vector<WordEntry>& evlist,
                      const Vector<WordEntry>& fvlist,
                      ostream& of2,
                      const Vector<WordIndex>& viterbi_alignment,
                      int pair_no, double alignment_score)
{
    WordIndex l, m;
    Vector<Vector<WordIndex> > translations(es.size()); 
    l = es.size() - 1;
    m = fs.size() - 1;
    
    Alignment a;
    for(WordIndex j = 1 ; j <= m ; j++)
        if(viterbi_alignment[j])
            a.insert(AlignmentPoint(viterbi_alignment[j]-1, j-1));
    
    ForceAligner::s_mutex.lock();
    if(ForceAligner::s_alignments.size() < pair_no)
        ForceAligner::s_alignments.resize(pair_no);
    ForceAligner::s_alignments[pair_no-1] = a;
    ForceAligner::s_mutex.unlock();    
}

/*********************************************************************************************/

GLOBAL_PARAMETER3(int,Model1_Iterations,"Model1_Iterations","NO. ITERATIONS MODEL 1","m1","number of iterations for Model 1",PARLEV_ITER,5);
GLOBAL_PARAMETER3(int,Model2_Iterations,"Model2_Iterations","NO. ITERATIONS MODEL 2","m2","number of iterations for Model 2",PARLEV_ITER,0);
GLOBAL_PARAMETER3(int,HMM_Iterations,"HMM_Iterations","mh","number of iterations for HMM alignment model","mh", PARLEV_ITER,5);
GLOBAL_PARAMETER3(int,Model3_Iterations,"Model3_Iterations","NO. ITERATIONS MODEL 3","m3","number of iterations for Model 3",PARLEV_ITER,5);
GLOBAL_PARAMETER3(int,Model4_Iterations,"Model4_Iterations","NO. ITERATIONS MODEL 4","m4","number of iterations for Model 4",PARLEV_ITER,5);
GLOBAL_PARAMETER3(int,Model5_Iterations,"Model5_Iterations","NO. ITERATIONS MODEL 5","m5","number of iterations for Model 5",PARLEV_ITER,0);
GLOBAL_PARAMETER3(int,Model6_Iterations,"Model6_Iterations","NO. ITERATIONS MODEL 6","m6","number of iterations for Model 6",PARLEV_ITER,0);
GLOBAL_PARAMETER(float, PROB_SMOOTH,"probSmooth","probability smoothing (floor) value ",PARLEV_OPTHEUR,1e-7);
GLOBAL_PARAMETER(float, MINCOUNTINCREASE,"minCountIncrease","minimal count increase",PARLEV_OPTHEUR,1e-7);
GLOBAL_PARAMETER2(int,Transfer_Dump_Freq,"TRANSFER DUMP FREQUENCY","t2to3","output: dump of transfer from Model 2 to 3",PARLEV_OUTPUT,0);
GLOBAL_PARAMETER2(bool,Verbose,"verbose","v","0: not verbose; 1: verbose",PARLEV_OUTPUT,0);
GLOBAL_PARAMETER(bool,Log,"log","0: no logfile; 1: logfile",PARLEV_OUTPUT,0);
GLOBAL_PARAMETER(double,P0,"p0","fixed value for parameter p_0 in IBM-3/4 (if negative then it is determined in training)",PARLEV_EM,-1.0);
GLOBAL_PARAMETER(double,M5P0,"m5p0","fixed value for parameter p_0 in IBM-5 (if negative then it is determined in training)",PARLEV_EM,-1.0);
GLOBAL_PARAMETER3(bool,Peg,"pegging","p","DO PEGGING? (Y/N)","0: no pegging; 1: do pegging",PARLEV_EM,0);
GLOBAL_PARAMETER(short,OldADBACKOFF,"adbackoff","",-1,0);
GLOBAL_PARAMETER2(unsigned int,MAX_SENTENCE_LENGTH,"ml","MAX SENTENCE LENGTH","maximum sentence length",0,MAX_SENTENCE_LENGTH_ALLOWED);
GLOBAL_PARAMETER(short, DeficientDistortionForEmptyWord,"DeficientDistortionForEmptyWord","0: IBM-3/IBM-4 as described in (Brown et al. 1993); 1: distortion model of empty word is deficient; 2: distoriton model of empty word is deficient (differently); setting this parameter also helps to avoid that during IBM-3 and IBM-4 training too many words are aligned with the empty word",PARLEV_MODELS,0);
GLOBAL_PARAMETER(int,restart,"restart","Restart training from a level,0: Normal restart, from model 1, 1: Model 1, 2: Model 2 Init (Using Model 1 model input and train model 2), 3: Model 2, (using model 2 input and train model 2), 4 : HMM Init (Using Model 1 model and train HMM), 5: HMM (Using Model 2 model and train HMM) 6 : HMM (Using HMM Model and train HMM), 7: Model 3 Init (Use HMM model and train model 3) 8: Model 3 Init (Use Model 2 model and train model 3) 9: Model 3, 10: Model 4 Init (Use Model 3 model and train Model 4) 11: Model 4 and on, ",PARLEV_INPUT,0);
GLOBAL_PARAMETER(bool,dumpCount,"dumpcount","Whether we are going to dump count (in addition to) final output?",PARLEV_OUTPUT,false);
GLOBAL_PARAMETER(bool,dumpCountUsingWordString,"dumpcountusingwordstring","In count table, should actual word appears or just the id? default is id",PARLEV_OUTPUT,false);
GLOBAL_PARAMETER(bool,ONLYALDUMPS,"ONLYALDUMPS","1: do not write any files",PARLEV_OUTPUT,0);
GLOBAL_PARAMETER(short,NCPUS,"NCPUS","Number of CPUS",PARLEV_EM,8);
GLOBAL_PARAMETER(short,CompactAlignmentFormat,"CompactAlignmentFormat","0: detailled alignment format, 1: compact alignment format ",PARLEV_OUTPUT,0);
GLOBAL_PARAMETER2(bool,NODUMPS,"NODUMPS","NO FILE DUMPS? (Y/N)","1: do not write any files",PARLEV_OUTPUT,0);
GLOBAL_PARAMETER(WordIndex, MAX_FERTILITY, "MAX_FERTILITY", "maximal fertility for fertility models", PARLEV_EM, 10);
short OutputInAachenFormat=0;
bool Transfer=TRANSFER;
bool Transfer2to3=0;
short NoEmptyWord=0;
bool FEWDUMPS=0;
Vector<map< pair<int,int>,char > > ReferenceAlignment;
bool useDict = false;
string Prefix, LogFilename, OPath, Usage, SourceVocabFilename, TargetVocabFilename;
string countPrefix;
Mutex logmsg_lock;
ofstream logmsg;
double LAMBDA=1.09;
string ReadTablePrefix;
void printGIZAPars(ostream&out) {}
double ErrorsInAlignment(const map< pair<int,int>,char >&reference,
		const Vector<WordIndex>&test, int l, int&missing, int&toomuch,
		int&eventsMissing, int&eventsToomuch, int pair_no) { return 1.0; }

