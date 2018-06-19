#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <functional>

#include "/opt/software/internal/yeastPolyA/progs23Nov/util.h"
#include "/opt/software/internal/yeastPolyA/progs23Nov/basicSeq.h"
#include "/opt/software/internal/yeastPolyA/progs23Nov/fastaSeq.h"

typedef std::map<std::string,std::string,std::less<std::string> > namedSeqMap;
// load fastaSeqs
// load bedfile
// extract and process all files
//
// Assume for now a small enough fastaSeq set that all can be loaded.  Presume existence of a samtools inddex file

int main(int argc,char *argv[]) {
  if(argc < 3) {
    std::cerr << "usage: " << argv[0] << " fastaFile bedFile [separator string]\n";
    exit(0);
  }
  std::string fafilename(argv[1]);
  std::string faifilename(argv[1]+std::string(".fai"));
  std::string bedfilename(argv[2]);
  std::string sep("_");
  if(argc>3) sep=argv[3];
  
  std::cerr << "Retrieving subsequences from file " << fafilename << " with index " << faifilename << " as defined in " << bedfilename << "\n";

  std::ifstream faifile(faifilename.c_str(),std::ios::in);
  std::ifstream fafile(fafilename.c_str(),std::ios::in);
  namedSeqMap faSeqs;
  while(!faifile.eof() && !fafile.eof()) {
    std::string seqID;
    int length,offset;
    char failine[1001];
    faifile.getline(failine,1000);
    strVec f=split()(std::string(failine));
    if(f.size()==5) {
      seqID=f[0];
      length=atoi(f[1].c_str());
      offset=atoi(f[2].c_str());
      //std::cerr << failine << " ==> " << seqID << ", " << length << ", " << offset << "\n";
      char faHeader[1001];
      fafile.getline(faHeader,1000);
      std::string seqID2=std::string(faHeader).substr(1,std::string(faHeader).length()-1);
      if(seqID != seqID2) { 
	std::cerr << "** Sequence ID mismatch: " << seqID << " != " << seqID2 << ". Exiting **\n";
	exit(1);
      }
      char *seqLine= new char[length+2];
      fafile.getline(seqLine,length+1);
      std::string seq(seqLine);
      if(length != seq.length()) {
	std::cerr << "** Sequence Length mismatch for sequence " << seqID << ": " << length << " != " << seq.length() << ". Exiting. **\n";
	exit(1);
      }
      faSeqs[seqID2]= seq;
      delete[] seqLine;
    }
  }
  std::cerr << "Loaded " << faSeqs.size() << " distinct sequences.\n";
  faifile.close();
  fafile.close();
  std::ifstream bedfile(bedfilename.c_str(),std::ios::in);
  int nOutOfRange=0;
  while(!bedfile.eof()) {
    char line[501];
    bedfile.getline(line,500);
    strVec f=split()(line);
    if(f.size()>=6) {
      std::string ID=f[0], strand=f[5], label=f[3];
      int start=atoi(f[1].c_str()),stop=atoi(f[2].c_str()), count=atoi(f[4].c_str());
      if(faSeqs.count(ID)) {
	if(stop>start && start <= faSeqs[ID].length() && stop<=faSeqs[ID].length()) {
	  std::string subseq= faSeqs[ID].substr(start,stop-start);
	  if(strand==std::string("-")) subseq= ReverseComplement()(subseq);
	  std::cout << label << '\t' << ID << "\t" << start << "\t" << stop << "\t" << strand << "\t" << count << "\t" << subseq
		    << '\t' << ID << sep << start << sep << stop << sep << strand << sep << count;
	  for(int i=6;i<f.size();++i) std::cout << '\t' << f[i];
	  std::cout << "\n";
	} else
	  ++nOutOfRange;
      }
    }
  }
  std::cerr << nOutOfRange << " out of range intervals\n"; 
  bedfile.close();
  return 0;
}
