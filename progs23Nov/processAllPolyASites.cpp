// program to convert polyA site counts to both site probability and also cumulative distributions
// written to tack on to the existing analysis, Joel Graber, December 2016
//
// We presume an input file of all sites from many processed polyA site sequences
// the file is tab-delimited with 11 or more columns, and one header row
//
// The first 10 columns of the input file are fixed with the following information:
//  1:#gene
//  2:chromo
//  3:position
//  4:strand
//  5:distToCDSstop
//  6:distToCDSstart
//  7:avgGcount
//  8:AGcount
//  9:Acount
//  10:Gcount
//
// All subsequent columns are the site counts for an individual sample, with a label in the first row of the file
//
// The input file is a required input parameter
// Additional parameters are
//   1: output file prefix: string prefix for output files, defaults to input file minus termnal ".txt"
//   2: ignore/scale/suppress multi-hit read:  0 for scale (default), 1 for ignore, 2 for suppress
//   3: AG count threshold: default to 8
//   4: CDS upstream cutoff distance:  default to 80
//   5: downstream cutoff distance:  default to 1000
// and for now MUST occur in this order
//
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include "/Users/knaggert/Desktop/pcf11_summer/progs23Nov/util.h"

int UScutoff=80,DScutoff=1000,AGthreshold=8,multiSwitch=0,paProbCutoff=30;
int ppCutoff=100;
double PseudoWeight=0.01;

typedef std::vector<double> dVec;
typedef std::vector<int> iVec;
typedef std::vector<std::string> strVec;
typedef std::map<std::string,int,std::less<std::string> > strIntMap;
typedef std::map<int,iVec> iIvecMap;
typedef std::map<int,strVec> iStrVecMap;
typedef std::map<std::string,iIvecMap> strIntIvecMap;
typedef std::map<std::string,iStrVecMap> strIntStrVecMap;

int main(int argc,char *argv[]) {
  if(argc<2) {
    std::cout << "usage: " << argv[0] << " polyAfilename [outputPrefix] [multiSwitch(0:1:2)=0] [AGcount=8] [UScutoff=80] [DScutoff=1000] [PseudoWeight=0.05] [paProbCutoff=30]\n";
    exit(0);
  }
  std::string ofPrefix(argv[1]);
  if(argc>2)
    ofPrefix=std::string(argv[2]);
  else { // remove ".txt" and any following text
    int p= ofPrefix.find(".txt");
    if(p>0) ofPrefix= ofPrefix.substr(0,p);
  }
  if(argc>3) { // must be 0, 1, or 2
    int tmp=atoi(argv[3]);
    if(tmp>=0 && tmp<3) multiSwitch=tmp;
  }
  if(argc>4) { // must be between 0 and 10 inclusive
    int tmp=atoi(argv[4]);
    if(tmp>=0 && tmp<11) AGthreshold=tmp;
  }
  if(argc>5) { // must be non-negative
    int tmp=atoi(argv[5]);
    if(tmp>=0) UScutoff=tmp;
  }
  if(argc>6) {// must be non-negative
    int tmp=atoi(argv[6]);
    if(tmp>=0) DScutoff=tmp;
  }
  if(argc>7) {// must be non-negative and <=1
    double tmp=atof(argv[7]);
    if(tmp>=0 && tmp<=1.0) PseudoWeight=tmp;
  }
  if(argc>8) {// must be non-negative
    int tmp=atoi(argv[8]);
    if(tmp>=0) paProbCutoff=tmp;
  }
  if(argc>9) {
    double tmp=atof(argv[9]);
    if(tmp>0.0) PseudoWeight=tmp;
  }
  std::cout << "#**************************************************************\n"
	    << "#* input file     : " << argv[1] << "\n"
	    << "#* output prefix  : " << ofPrefix << "\n"
	    << "#* multiSwitch    : " << multiSwitch << "\n"
	    << "#* AGthreshold    : " << AGthreshold << "\n"
	    << "#* US cutoff      : " << UScutoff << "\n"
	    << "#* DS cutoff      : " << DScutoff << "\n"
	    << "#* pA prob cutoff : " << paProbCutoff << "\n"
	    << "#* PseudoWeight   : " << PseudoWeight << "\n"
	    << "#**************************************************************\n";

  std::ifstream infile(argv[1],std::ios::in);
  char lineIn[10001];
  infile.getline(lineIn,10000);
  strVec headers=split()(std::string(lineIn));
  if(headers.size()<11) {
    std::cerr << "Insufficient columns in input file: " << argv[1] << ", must be >10, found" << headers.size() << ", exiting\n";
    exit(0);
  }
  strVec sampleLabels; for(int i=10;i<headers.size();++i) sampleLabels.push_back(headers[i]);
  int Nsamples=sampleLabels.size();
  int lineCtr=1,validDataLines=0,keptDataLines=0,cut5lines=0,cut3lines=0,cutAGlines=0,cutMultilines=0;
  strIntIvecMap GeneCountMap;
  strIntStrVecMap GeneMetadataMap;
  while(!infile.eof()) {
    infile.getline(lineIn,10000);
    strVec f=split()(std::string(lineIn));
    if(f.size()>0) {
      if(f.size()!=headers.size()) {
	std::cerr << "Column count mismatch at line " << ++lineCtr << ", expected " << headers.size() << ", read " << f.size() << "; exiting.";
	exit(0);
      }
      ++validDataLines;
      // check for AGthreshold and distance cutoffs
      bool AGokay= (atoi(f[7].c_str())<=AGthreshold);
      bool cut5prOkay= (atoi(f[5].c_str()) >= (-1*UScutoff));
      bool cut3prOkay= (atoi(f[4].c_str()) <= DScutoff);
      bool cutMultiOkay= (multiSwitch<2 || atof(f[6].c_str())<1.5);
      if(AGokay && cut5prOkay && cut3prOkay && cutMultiOkay) {
	++keptDataLines;
	// line is okay, so now process it.  The distance to CDS stop is the index and will properly sort the sites for both W and C genes
	int DistToCDSstop= atoi(f[4].c_str());
	strVec metadata; for(int i=0;i<10;++i) metadata.push_back(f[i]);
	iVec countdata; for(int i=10;i<f.size();++i) countdata.push_back(atoi(f[i].c_str()));
	if(!GeneCountMap.count(f[0])) {
	  GeneCountMap[f[0]]= iIvecMap();
	  GeneMetadataMap[f[0]]= iStrVecMap();
	}
	GeneCountMap[f[0]][DistToCDSstop]= countdata;
	GeneMetadataMap[f[0]][DistToCDSstop]= metadata;
      } else {
	if(!AGokay) ++cutAGlines;
	if(!cut5prOkay) ++cut5lines;
	if(!cut3prOkay) ++cut3lines;
	if(!cutMultiOkay) ++cutMultilines;
      }
    }
    ++lineCtr;
  }
  infile.close();
  std::cout << "# samples found   : " << Nsamples << "\n"
	    << "# total lines     : " << lineCtr << "\n"
	    << "# valid lines     : " << validDataLines << "\n"
	    << "# kept lines      : " << keptDataLines << "\n"
	    << "# AG-fail lines   : " << cutAGlines << "\n"
	    << "# 5'-cut lines    : " << cut5lines << "\n"
	    << "# 3'-cut lines    : " << cut3lines << "\n"
	    << "# Multi-cut lines :" << cutMultilines << "\n"
	    << "# Genes counted: " << GeneMetadataMap.size() << "\n";

  std::string cFilename(ofPrefix+std::string("_cumPa.txt")),pFilename(ofPrefix+std::string("_paProb.txt")),gFilename(ofPrefix+std::string("_paGene.txt")),
    gcFilename(ofPrefix+std::string("_geneCounts.txt")),flcFilename(ofPrefix+std::string("_flGeneCounts.txt"));;
  std::ofstream cfile(cFilename.c_str(),std::ios::out);
  std::ofstream pfile(pFilename.c_str(),std::ios::out);
  std::ofstream gfile(gFilename.c_str(),std::ios::out);
  std::ofstream gcfile(gcFilename.c_str(),std::ios::out);
  std::ofstream flcfile(flcFilename.c_str(),std::ios::out);
  
  strIntIvecMap::iterator cIt=GeneCountMap.begin();
  strIntStrVecMap::iterator mdIt=GeneMetadataMap.begin();
  for(int i=0;i<headers.size();++i) {
    if(i) { cfile << '\t'; pfile << '\t'; }
    cfile << headers[i];
    pfile << headers[i];
  }
  cfile << "\tall\n";
  pfile << "\tall\n";
  // headers for the gene file are different- just the gene name (header[0]), then the sample id three times over with suffixes
  gfile << headers[0] << "\tSampleTotal\tUTR3count\tppCount";
  for(int i=10;i<headers.size();++i) gfile << '\t' << headers[i] << "_AvgL3pUTR";
  for(int i=10;i<headers.size();++i) gfile << '\t' << headers[i] << "_fracTrunc";
  for(int i=10;i<headers.size();++i) gfile << '\t' << headers[i] << "_fracProPrx";
  gfile << "\n";
  // headers for the gene counts files
  gcfile << headers[0];
  for(int i=10;i<headers.size();++i) gcfile << '\t' << headers[i];
  gcfile << "\n";
  flcfile << headers[0];
  for(int i=10;i<headers.size();++i) flcfile << '\t' << headers[i];
  flcfile << "\n";
  //
  typedef std::map <int,double> idMap;
  // loop through all genes; both GeneCountMap and GeneMetadataMap are indexed by a string reprensting the gene name
  //   meaning that (*cIt).first and (*mdIt).first are the same and are both the gene name.
  while(cIt != GeneCountMap.end() && mdIt != GeneMetadataMap.end()) {
    //first get the sample totals for all, including scaling if multiSwitch==0
    dVec sampleTotals(Nsamples+1,0.0);
    //vectors to hold gene level data:
    //    UTR3count is the number of sites past the stop codon
    //    UTR3sum is the weighted sum of sites past the stop codon
    //    ppCount is the number of sites within threshold (100nt by default) of the start codon
    //  the truncated count can be calculated as the complement of the UTR3count
    dVec UTR3count(Nsamples+1,0.0),UTR3sum(Nsamples+1,0.0),ppCount(Nsamples+1,0.0);
    idMap allSampleCount;
    //loop through each position (the int that is the key of the map that is iterated through by it)
    for(iIvecMap::iterator it=(*cIt).second.begin();it!=(*cIt).second.end();++it) {
      double scale=1.0;
      if(multiSwitch==0) scale= atof((*mdIt).second[(*it).first][6].c_str());
      allSampleCount[(*it).first]=0.0;
      for(int i=0;i<(*it).second.size();++i) {
	double val=static_cast<double>((*it).second[i])/scale;
	sampleTotals[i] += val;
	sampleTotals[Nsamples] += val;
	allSampleCount[(*it).first] += val;
	// collecting gene-level info
	if((*it).first>=0) { // distToCDSstop >=0 means site is in 3'-UTR
	  UTR3count[i] += val;
	  UTR3count[Nsamples] += val;
	  UTR3sum[i] += val*(*it).first;
	  UTR3sum[Nsamples] += val*(*it).first;
	} else if(atoi((*mdIt).second[(*it).first][5].c_str())<=ppCutoff) { // within ppCutoff of startCodon 
	  ppCount[i] += val;
	  ppCount[Nsamples] += val;
	}
      }
    }
    // write the gene data out to gfile
    gfile << (*cIt).first << '\t' << sampleTotals[Nsamples] << '\t' << UTR3count[Nsamples] << '\t' << ppCount[Nsamples];
    // average 3'-UTR length (counting only those sites that land in the 3'-UTR
    for(int i=0;i<Nsamples;++i)
      if((UTR3count[i]+PseudoWeight*UTR3count[Nsamples])>0.0)
	gfile << '\t' << (UTR3sum[i]+PseudoWeight*UTR3sum[Nsamples])/(UTR3count[i]+PseudoWeight*UTR3count[Nsamples]);
      else
	gfile << '\t' << "-1";
    // fraction of transcripts that are truncated
    for(int i=0;i<Nsamples;++i) gfile << '\t' << (1.0 - (UTR3count[i]+PseudoWeight*UTR3count[Nsamples])/(sampleTotals[i]+PseudoWeight*sampleTotals[Nsamples]));
    // fraction of transcripts that are within ppCutoff of startCodon
    for(int i=0;i<Nsamples;++i) gfile << '\t' << ((ppCount[i]+PseudoWeight*ppCount[Nsamples])/(sampleTotals[i]+PseudoWeight*sampleTotals[Nsamples]));    
    gfile << '\n';
    // gene count files
    gcfile << (*cIt).first; flcfile << (*cIt).first;
    for(int i=0;i<Nsamples;++i) gcfile << '\t' << sampleTotals[i];
    for(int i=0;i<Nsamples;++i) flcfile << '\t' << UTR3count[i];
    gcfile << "\n"; flcfile << "\n";
    //
    // *** back to site by site calculation ***
    // now we have all (scaled) total counts for each sample for this gene, and also all summed counts at each position
    // first make the cumulative counts
    dVec runningCounts(Nsamples+1,0.0);
    for(iIvecMap::iterator it=(*cIt).second.begin();it!=(*cIt).second.end();++it) {
      // first calculate the remaining counts for each sample
      dVec paProb(Nsamples+1,0.0),remCounts(Nsamples+1,0.0);
      for(int i=0;i<=Nsamples;++i) remCounts[i]=sampleTotals[i]-runningCounts[i];
      // now get the scaled counts for this site
      double scale=1.0,siteSum=0.0;
      if(multiSwitch==0) scale= atof((*mdIt).second[(*it).first][6].c_str());
      for(int i=0;i<(*it).second.size();++i) {
	double val=static_cast<double>((*it).second[i])/scale;
	runningCounts[Nsamples] += val;
	runningCounts[i] += val;
	siteSum+=val;
      }	
      for(int i=0;i<(*it).second.size();++i) {
	double val=static_cast<double>((*it).second[i])/scale;
	// 05Jan2017- switching to using total counts across all samples as a weighted basis for a Beta Prior
	//  as as consequence, probabilities are calculated for any site in which the sum of remaining counts across
	//  all samples exceeds the threshold
	//  old code: if(remCounts[i]>=paProbCutoff) paProb[i]= val/remCounts[i];
	if(remCounts[Nsamples]>=paProbCutoff)
	  paProb[i]= (val+PseudoWeight*(siteSum+1.0))/(remCounts[i]+PseudoWeight*(remCounts[Nsamples]+1));
      }
      if(remCounts[Nsamples]>=paProbCutoff) paProb[Nsamples]=siteSum/remCounts[Nsamples];
      if(siteSum > 0.0) {
	// write the metadata out to both files
	for(int i=0;i<10;++i) { // metadata first
	  if(i) { cfile << '\t'; pfile << '\t'; }
	  cfile << (*mdIt).second[(*it).first][i];
	  pfile << (*mdIt).second[(*it).first][i];
	}
	// cumulative counts
	for(int i=0;i<=Nsamples;++i)
	  cfile << '\t' << (runningCounts[i]/sampleTotals[i]);
	// probability at this site
	for(int i=0;i<=Nsamples;++i)
	  pfile << '\t' << paProb[i];
	cfile << '\n';
	pfile << '\n';
      }
    }
    ++cIt;  ++mdIt;
  }
  cfile.close();
  pfile.close();
  gfile.close();
  gcfile.close();
  flcfile.close();
  return 0;
}
    
    
