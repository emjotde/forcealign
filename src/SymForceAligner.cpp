#include "SymForceAligner.h"

SymForceAligner::SymForceAligner(std::string src, std::string trg, std::string path, size_t CPUS)
: m_mode(GrowDiagFinalAnd), m_diagonal(true), m_final(true), m_bothuncovered(true),
  m_srcTrgAligner(src,trg,path, CPUS), m_trgSrcAligner(trg, src, path, CPUS)
{ }

void SymForceAligner::setMode(Mode mode) {
    m_mode = mode;
}

Alignment SymForceAligner::cUnion(int *a, int m, int* b, int n){
    Alignment out;
    for (int j=1;j<=m;j++)
        if (a[j]) 
            out.insert(AlignmentPoint(a[j]-1, j-1));

    for (int i=1;i<=n;i++)
        if (b[i] && a[b[i]]!=i)
            out.insert(AlignmentPoint(i-1, b[i]-1));
              
    return out;
}

Alignment SymForceAligner::cIntersection(int *a, int m, int* b, int n){
    Alignment out;
    for (int j = 1; j <= m; j++)
        if (a[j] && b[a[j]] == j)
            out.insert(AlignmentPoint(a[j]-1, j-1));
    return out;
}

Alignment SymForceAligner::cGrow(int *a, int m, int *b, int n) {
    int* fa = new int[MAX_M + 1];
    int* ea = new int[MAX_N + 1];
    
    int** A = new int*[MAX_N + 1];
    for(int i=1; i <= MAX_N; i++)
        A[i]=new int[MAX_M + 1];
    
    std::vector<std::pair <int,int> > neighbors;
    std::pair<int,int> entry;
 
    neighbors.push_back(make_pair(-1,-0));
    neighbors.push_back(make_pair(0,-1));
    neighbors.push_back(make_pair(1,0));
    neighbors.push_back(make_pair(0,1));
 
    if (m_diagonal){
        neighbors.push_back(make_pair(-1,-1));
        neighbors.push_back(make_pair(-1,1));
        neighbors.push_back(make_pair(1,-1));
        neighbors.push_back(make_pair(1,1));    
    }

    int i,j,o;
    //covered foreign and english positions      
    memset(fa, 0 ,(m + 1)*sizeof(int));
    memset(ea, 0 ,(n + 1)*sizeof(int));
 
    //matrix to quickly check if one point is in the symmetric
    //alignment (value=2), direct alignment (=1) and inverse alignment   

    for (int i = 1; i <= n; i++)
        memset(A[i], 0, (m + 1) * sizeof(int));
 
    std::set<std::pair<int,int> > currentpoints; //symmetric alignment
    std::set<std::pair<int,int> > unionalignment; //union alignment
 
    std::pair<int,int> point; //variable to store points
    std::set<std::pair<int,int> >::const_iterator k; //iterator over sets
 
    //fill in the alignments
    for(j=1; j<=m; j++){
       if(a[j]){
          unionalignment.insert(make_pair(a[j], j));
          if (b[a[j]] == j){ 
             fa[j] = 1;
             ea[a[j]] = 1;
             A[a[j]][j] = 2;   
             currentpoints.insert(make_pair(a[j], j));
          }         
          else 
             A[a[j]][j] = -1;
       }
    }
    
    for(i=1; i<=n; i++) 
       if (b[i] && a[b[i]] !=i ){ //not intersection
          unionalignment.insert(make_pair(i, b[i])); 
          A[i][b[i]] = 1;
       } 

    int added=1;
    while(added) {
       added = 0;
       ///scan the current alignment
       for (k=currentpoints.begin();k!=currentpoints.end();k++){
          //cout << "{"<< (k->second)-1 << "-" << (k->first)-1 << "}";
          for (o=0;o<neighbors.size();o++){
             //cout << "go over check all neighbors\n";
             point.first=k->first+neighbors[o].first;      
             point.second=k->second+neighbors[o].second; 
             //cout << point.second-1 << " " << point.first-1 << "\n";
             //check if neighbor is inside 'matrix'
             if (point.first>0 && point.first <=n && point.second>0 && point.second<=m)
                //check if neighbor is in the unionalignment alignment
                if (b[point.first]==point.second || a[point.second]==point.first){
                   //cout << "In unionalignment ";cout.flush();
                   //check if it connects at least one uncovered word
                   if (!(ea[point.first] && fa[point.second]))
                   {
                      //insert point in currentpoints!
                      currentpoints.insert(point);
                      A[point.first][point.second]=2;
                      ea[point.first]=1; fa[point.second]=1;
                      added=1;
                      //cout << "added grow: " << point.second-1 << "-" << point.first-1 << "\n";cout.flush();
                   }
                }
          }
       }
    }
       
    if (m_final){
       for (k=unionalignment.begin();k!=unionalignment.end();k++)
          if (A[k->first][k->second]==1)
          {            
             point.first=k->first;point.second=k->second;
             //one of the two words is not covered yet
             //cout << "{" << point.second-1 << "-" << point.first-1 << "} ";
             if ((m_bothuncovered &&  !ea[point.first] && !fa[point.second]) ||
                (!m_bothuncovered && !(ea[point.first] && fa[point.second])))
             {
                //add it!
                currentpoints.insert(point);
                A[point.first][point.second]=2;
                //keep track of new covered positions                
                ea[point.first]=1;fa[point.second]=1;
                
                //added=1;
                //cout << "added final: " << point.second-1 << "-" << point.first-1 << "\n";
             }
          }
             
             for (k=unionalignment.begin();k!=unionalignment.end();k++)
                if (A[k->first][k->second]==-1)
                {            
                   point.first=k->first;point.second=k->second;
                   //one of the two words is not covered yet
                   //cout << "{" << point.second-1 << "-" << point.first-1 << "} ";
                   if ((m_bothuncovered &&  !ea[point.first] && !fa[point.second]) ||
                       (!m_bothuncovered && !(ea[point.first] && fa[point.second])))
                   {
                      //add it!
                      currentpoints.insert(point);
                      A[point.first][point.second]=2;
                      //keep track of new covered positions                
                      ea[point.first]=1;fa[point.second]=1;
                      
                      //added=1;
                      //cout << "added final: " << point.second-1 << "-" << point.first-1 << "\n";
                   }
                }
    }
    
    Alignment out;
    for (k=currentpoints.begin();k!=currentpoints.end();k++) {
        out.insert(AlignmentPoint(k->first-1, k->second-1));
    }

    delete[] A;
    delete[] ea;
    delete[] fa;
    
    return out;
}

Alignment invert(Alignment a) {
    Alignment b;
    for(Alignment::iterator it = a.begin(); it != a.end(); it++)
        b.insert(AlignmentPoint(it->second, it->first));
    return b;
} 

Alignment SymForceAligner::symmetrize(Alignment& a_, Alignment& b_) {
    int a[MAX_M], b[MAX_N], m, n;
    std::fill_n(a, MAX_M, 0);
    std::fill_n(b, MAX_N, 0);
    
    m = 0;
    for(Alignment::iterator it = a_.begin(); it != a_.end(); it++) {
        a[it->second+1] = it->first+1;
        if(it->second+1 > m)
            m = it->second+1;
    }
    
    n = 0;
    for(Alignment::iterator it = b_.begin(); it != b_.end(); it++) {
        b[it->second+1] = it->first+1;
        if(it->second+1 > n)
            n = it->second+1;
    }
    
    switch(m_mode) {
      case Src2Trg:
        return a_;
      case Trg2Src:
        return invert(b_);
      case Intersection:
        return cIntersection(a, m, b, n);
      case Union:
        return cUnion(a, m, b, n);
      case Grow:
        m_diagonal = false; m_final = false; m_bothuncovered = false;
        return cGrow(a, m, b, n);
      case GrowDiag:
        m_diagonal = true; m_final = false; m_bothuncovered = false;
        return cGrow(a, m, b, n);
      case GrowDiagFinal:
        m_diagonal = true; m_final = true; m_bothuncovered = false;
        return cGrow(a, m, b, n);
      case GrowDiagFinalAnd:
        m_diagonal = true; m_final = true; m_bothuncovered = true;
        return cGrow(a, m, b, n);
    }
    
    return cGrow(a, m, b, n);
}

void SymForceAligner::addSentence(std::string e, std::string f) {
    m_srcTrgAligner.addSentence(e, f);
    m_trgSrcAligner.addSentence(f, e);
}

void SymForceAligner::alignCorpus() {
    m_alignments.clear();
  
    m_srcTrgAligner.alignCorpus();
    m_trgSrcAligner.alignCorpus();
    
    Alignments::iterator it1 = m_srcTrgAligner.getAlignments().begin();
    Alignments::iterator it2 = m_trgSrcAligner.getAlignments().begin();
    
    while(it1 != m_srcTrgAligner.getAlignments().end() && it2 != m_trgSrcAligner.getAlignments().end()) {
        Alignment a1 = *it1;
        Alignment a2 = *it2;
        
        Alignment a = symmetrize(a1, a2);
        m_alignments.push_back(a);
        it1++; it2++;
    }    
}

Alignment& SymForceAligner::alignSentence(std::string e, std::string f) {  
    clearAll();
    addSentence(e, f);
    alignCorpus();
    return m_alignments[0];
}

std::string SymForceAligner::alignSentenceStr(std::string e, std::string f) {
    Alignment a = alignSentence(e, f);

    std::stringstream ss;
    for(Alignment::iterator it = a.begin(); it != a.end(); it++) {
        if(it != a.begin())
            ss << " ";
        ss << it->first << "-" << it->second;
    }
    return ss.str();
}

void SymForceAligner::clearAll() {
    m_alignments.clear();
    m_srcTrgAligner.clearAll();
    m_trgSrcAligner.clearAll();
}

Alignments& SymForceAligner::getAlignments() { return m_alignments; }

std::string SymForceAligner::getAlignmentsStr() {
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

