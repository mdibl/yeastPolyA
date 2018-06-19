#include <iostream>
#include <iomanip>
#include <fstream>
#include "/Users/jhgraber/prog/cpp_core/util.h"

#include <string>
#include <map>
#include <functional>
#include <algorithm>
#include <vector>


//typedef std::vector<double> dVec;
typedef std::map<std::string,dVec> StrToDVecMap;

// vector of counts to be added for each entity (either in the gene or in the "feature" category)
//  0  total counts
//  1  non-Arich counts
//  2  5' counts
//  3  internal counts
//  4  3' counts
//  5  non-Arich 5' counts
//  6  non-Arich internal counts
//  7  non-Arich 3' counts
//  8  Scaled total counts (dividing counts by genomic location counts)
//  9  Scaled non-Arich counts
// 10  Scaled 5' counts
// 11  Scaled internal counts
// 12  Scaled 3' counts
// 13  Scaled non-Arich 5' counts
// 14  Scaled non-Arich internal counts
// 15  Scaled non-Arich 3' counts

// input parameter- number of As and Gs in the DS 10nt sequence that triggers "Arich designation" -- default to 8
int ArichTH=8;
const dVec zeroVec(16,0.0);

// Presumed input line configuation:
//  0: sample label             [UA00e1]
//  1: chromosome               [chrIII]
//  2: site                     [100335]
//  3: strand                   [+]
//  4: count                    [3]
//  5: seqCount                 [3]
//  6: genAlnCount              [1]
//  7: upstream seq             [ATAATAATCCATCCAATTCTATTCC]
//  8: downstream seq           [AAAAATCGAGAAGCCCCCAGCATTC]
//  9: ds AG count              [8]
// 10: ds A count               [6]
// 11: ds G count               [2]
// 12: seqlabel                 [chrIII_100310_100360_+_3]
// 13: gene-site loc            [int]
// 14: dist to closest gene end [856]
// 15: closest gene             [YCL014W]
// 16: closest gene start       [96281]
// 17: closest gene end         [101191]
// 18: closest gene type        [gene]
// 19: feature-site loc         [f_int]
// 20: dist to closest feat end [856]
// 21: closest feature          [YCL014W]
// 22: closest feature start    [96281]
// 23: closest feature end      [101191]
// 24: closest feature type     [gene]
 
int main(int argc, char *argv[]) {
  if(argc <3) {
    std::cout << "Usage: " << argv[0] << " inputFile outputPrefix [ArichTH]\n";
    exit(0);
  }
  if(argc > 3) {
    int th=atoi(argv[3]);
    if(th>=0 && th<=10) ArichTH=th;
  }

  // Load all gene and feature labels to make uniform arrays across all files
  std::ifstream gListFile("/Users/jhgraber/data/yeast/saccer3_simple_features_allGenes.txt",std::ios::in);
  strVec geneList;
  while(!gListFile.eof()) {
    char line[1001];
    gListFile.getline(line,1000);
    if(std::string(line).length()>0)
      geneList.push_back(std::string(line));
  }
  gListFile.close();
  std::cerr << "Loaded " << geneList.size() << " gene names\n";
  std::ifstream fListFile("/Users/jhgraber/data/yeast/saccer3_simple_features_allFeatures.txt",std::ios::in);
  strVec featList;
  while(!fListFile.eof()) {
    char line[1001];
    fListFile.getline(line,1000);
    if(std::string(line).length()>0)
      featList.push_back(std::string(line));
  }
  fListFile.close();
  std::cerr << "Loaded " << featList.size() << " feature names\n";
  // load the count maps
  StrToDVecMap GeneCounts,FeatureCounts;
  for(int i=0;i<geneList.size();++i) GeneCounts[geneList[i]]= zeroVec;
  for(int i=0;i<featList.size();++i) FeatureCounts[featList[i]]= zeroVec;
  std::cerr << "GeneCounts= " << GeneCounts.size() << ", FeatCounts=" << FeatureCounts.size() << '\n';
  // process the input file
  std::ifstream infile(argv[1],std::ios::in);
  int ctr=0;
  while(!infile.eof()) {
    if(!(++ctr % 10000)) { std::cout << "."; std::cout.flush(); }
    char lineIn[5001];
    infile.getline(lineIn,5000);
    strVec f=split()(std::string(lineIn));
    if(f.size()>=25) {
      double totCount= atof(f[4].c_str()), scaledCount= atof(f[4].c_str())/atof(f[6].c_str());
      int agCount= atoi(f[9].c_str());
      std::string gID=f[15],fID=f[21];
      // first process the gene data
      GeneCounts[gID][0]    += totCount; GeneCounts[gID][8]    += scaledCount;
      FeatureCounts[fID][0] += totCount; FeatureCounts[fID][8] += scaledCount;
      // gene by location
      if(f[13] == std::string("ext_5"))    {
	GeneCounts[gID][2] += totCount; GeneCounts[gID][10] += scaledCount;
      } else if(f[13] == std::string("int")) {
	GeneCounts[gID][3] += totCount; GeneCounts[gID][11] += scaledCount;
      } else {
	GeneCounts[gID][4] += totCount; GeneCounts[gID][12] += scaledCount;
      }
      // feature by location
      if(f[19] == std::string("f_ext_5"))    {
	FeatureCounts[fID][2] += totCount; FeatureCounts[fID][10] += scaledCount;
      } else if(f[19] == std::string("f_int")) {
	FeatureCounts[fID][3] += totCount; FeatureCounts[fID][11] += scaledCount;
      } else {
	FeatureCounts[fID][4] += totCount; FeatureCounts[fID][12] += scaledCount;
      }
      //
      // now process again if the ds sequence is not AG-rich
      if(agCount < ArichTH) {
	GeneCounts[gID][1]    += totCount; GeneCounts[gID][9]    += scaledCount;
	FeatureCounts[fID][1] += totCount; FeatureCounts[fID][9] += scaledCount;
	// gene by location                                       	
	if(f[13] == std::string("ext_5"))    {
	  GeneCounts[gID][5] += totCount; GeneCounts[gID][13] += scaledCount;
	} else if(f[13] == std::string("int")) {
	  GeneCounts[gID][6] += totCount; GeneCounts[gID][14] += scaledCount;
	} else {
	  GeneCounts[gID][7] += totCount; GeneCounts[gID][15] += scaledCount;
	}
	// feature by location
	if(f[19] == std::string("f_ext_5"))    {
	  FeatureCounts[fID][5] += totCount; FeatureCounts[fID][13] += scaledCount;
	} else if(f[19] == std::string("f_int")) {
	  FeatureCounts[fID][6] += totCount; FeatureCounts[fID][14] += scaledCount;
	} else {
	  FeatureCounts[fID][7] += totCount; FeatureCounts[fID][15] += scaledCount;
	}
      }
    }
  } // while not end of file
  infile.close(); 
  int geneCount=GeneCounts.size(),featCount=FeatureCounts.size();
  int g3q= 3*geneCount/4, f3q=3*featCount/4;
  std::cerr << "\n* Gene Map Size= " << geneCount << "\n* Feat Map Size= " << featCount << "\n";
  std::cerr << "\n* Gene 3rd Quar= " << g3q << "\n* Feat 3rd Quar= " << f3q << "\n";
  std::cerr << "Calculating 3rd Quartiles..."; std::cerr.flush();
  // calculate third quartiles
  dVec g3rdQ=zeroVec,f3rdQ=zeroVec;
  for(int i=0;i<g3rdQ.size();++i) {
    //std::cerr << " " << i; std::cerr.flush();
    dVec gtemp;
    for(StrToDVecMap::iterator it=GeneCounts.begin();it!=GeneCounts.end();++it) if((*it).second[i]>0) gtemp.push_back((*it).second[i]); 
    //std::cerr << gtemp.size() << 'a'; std::cerr.flush();
    if(gtemp.size()) std::sort(gtemp.begin(),gtemp.end());
    g3rdQ[i]= gtemp[3*gtemp.size()/4];
    //std::cerr << 'a'; std::cerr.flush();
    //
    dVec ftemp;
    for(StrToDVecMap::iterator it=FeatureCounts.begin();it!=FeatureCounts.end();++it) if((*it).second[i]>0) ftemp.push_back((*it).second[i]);
    //std::cerr << 'b'; std::cerr.flush();
    if(ftemp.size()) std::sort(ftemp.begin(),ftemp.end());
    f3rdQ[i]= ftemp[3*ftemp.size()/4];
    //std::cerr << 'b'; std::cerr.flush();
  }

  // write output files
  std::string gfname(std::string(argv[2])+std::string(".gene.txt"));
  std::cerr << "done.\n\n Writing Gene Count File " << gfname << "..."; std::cerr.flush();
  std::ofstream gfile(gfname.c_str(),std::ios::out);
  gfile << "#gene\ttotal\tnonA\t5p\tint\t3p\tnonA5p\tnonAint\tnonA3p\tscaled\tsc_nonA\tsc_5p\tsc_int\tsc_3p\tsc_nonA5p\tsc_nonAint\tsc_nonA3p\n";
  for(StrToDVecMap::iterator it=GeneCounts.begin();it!=GeneCounts.end();++it) {
    gfile << (*it).first;
    for(int i=0;i<(*it).second.size();++i) gfile << '\t' << (*it).second[i];
    gfile << '\n';
  }
  gfile << "Third_Quartile";
  for(int i=0;i<g3rdQ.size();++i) gfile << '\t' << g3rdQ[i];
  gfile << '\n';
  gfile.close();
  std::string ffname(std::string(argv[2])+std::string(".feature.txt"));
  std::cerr << "\n Writing Feature Count File " << ffname << "..."; std::cerr.flush();
  std::ofstream ffile(ffname.c_str(),std::ios::out);
  ffile << "#gene\ttotal\tnonA\t5p\tint\t3p\tnonA5p\tnonAint\tnonA3p\tscaled\tsc_nonA\tsc_5p\tsc_int\tsc_3p\tsc_nonA5p\tsc_nonAint\tsc_nonA3p\n";
  for(StrToDVecMap::iterator it=FeatureCounts.begin();it!=FeatureCounts.end();++it) {
    ffile << (*it).first;
    for(int i=0;i<(*it).second.size();++i) ffile << '\t' << (*it).second[i];
    ffile << '\n';
  }    
  ffile << "Third_Quartile";
  for(int i=0;i<f3rdQ.size();++i) ffile << '\t' << f3rdQ[i];
  ffile << '\n';
  ffile.close();
  std::cerr << "done.\n";

  return 0;
}
// vector of counts to be added for each entity (either in the gene or in the "feature" category)
//  0  Total counts
//  1  non-Arich counts    
//  2  5' counts    
//  3  internal counts    
//  4  3' counts      
//  5  non-Arich 5' counts    
//  6  non-Arich internal counts    
//  7  non-Arich 3' counts          
//  8  Scaled total counts (dividing counts by genomic location counts)    
//  9  Scaled non-Arich counts     
// 10  Scaled 5' counts         
// 11  Scaled internal counts       
// 12  Scaled 3' counts                
// 13  Scaled non-Arich 5' counts              
// 14  Scaled non-Arich internal counts           
// 15  Scaled non-Arich 3' counts   
