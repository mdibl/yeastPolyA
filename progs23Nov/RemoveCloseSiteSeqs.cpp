// assume that we have a fasta file that has labels of the form >chr_start_stop_str_pos, followed by a single line of sequence
//  we just need to get rid of any that are two close together to not skew downstream sequence analysis
//  take a parameter for distance separated- chromosome and strant must also match.  The first one found will be taken.
// for now, no error checking on the input file, assume it will be correct (01 Feb 2017)
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

int thDist=20;  // the closest spacing allowed.

typedef std::vector<int> iVec;
typedef std::vector<std::string> strVec;

int main(int argc, char *argv[]) {

  if(argc<2) {
    std::cout << "usage: " << argv[0] << " filename [int:distance]\n";
    exit(0);
  }
  if(argc > 2) {
    int t= atoi(argv[2]);
    if(t>0) thDist=t;
  }

  std::ifstream infile(argv[1],std::ios::in);
  iVec posVec;
  strVec chrStrVec,headVec,seqVec;
  std::string infilename(argv[1]);
  std::string outfilename(infilename.substr(0,infilename.length()-3)+std::string("_filtered.fa"));
  std::cerr << "Processing fasta file: " << argv[1] << " with thDist= " << thDist << ", writing output to: " << outfilename << "\n";

  int c=0,f=0;
  while(!infile.eof()) {
    char lineIn[10001];
    infile.getline(lineIn,10000);
    std::string head(lineIn);
    if(head.length()>0) {
      int p1= head.find("_"), p2=head.rfind("_");
      int pos= atoi(head.substr(p2+1,head.length()).c_str());
      std::string chrStr= head.substr(0,p1)+head.substr(p2-1,1);
      //std::cout << head << " ==> " << chrStr << "\t" << pos << "\n";
      infile.getline(lineIn,10000);
      std::string seq(lineIn);
      bool filter=false;
      for(int i=0;i<chrStrVec.size()&& !filter; ++i)
	filter = (chrStr==chrStrVec[i] && abs(pos-posVec[i]) <=thDist);
      if(!filter) {
	headVec.push_back(head);
	seqVec.push_back(seq);
	posVec.push_back(pos);
	chrStrVec.push_back(chrStr);
	++c;
      } else
	++f;
    }
  }
  infile.close();
  std::cerr << "kept " << c << ", filtered " << f << "\n";
  std::ofstream outfile(outfilename.c_str(),std::ios::out);
  for(int i=0;i<headVec.size();++i)
    outfile << headVec[i] << '\n' << seqVec[i] << '\n';
  outfile.close();

  return 0;
}
      
      


    
	
