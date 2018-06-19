#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <algorithm>

#include "/opt/software/internal/yeastPolyA/progs23Nov/util.h"

typedef std::vector<int> iVec;
// string to integer map, for sample to counts- at the end, all entries will be filled including 0 entries
typedef std::map<std::string,int> strIntMap;
// integer (for position) to a map of sample-to-count (string to int)
typedef std::map<int,strIntMap> intStrIntMap;
// also need a structure to hold info about the site (multi-alignment and AGcounts)
class siteInfo {
public:
  int agCounts[3];
  double gCount;
  siteInfo():agCounts(),gCount() {}
  siteInfo(const siteInfo& o) { gCount=o.gCount; for(int i=0;i<3;++i) agCounts[i]=o.agCounts[i]; }
  siteInfo(int *agc,double g) { gCount=g; for(int i=0;i<3;++i) agCounts[i]=agc[i]; }
  siteInfo& operator=(const siteInfo& o) { gCount=o.gCount; for(int i=0;i<3;++i) agCounts[i]=o.agCounts[i]; return *this; }
};
// map of int (position) to siteInfo structure
typedef std::map<int,siteInfo> intSiteInfoMap;
// struct to hold everything by gene- // ID not needed because it will be the key for the map that holds it
class geneData {
public:
  std::string chromo,strand; // chromosome and strand
  int start,stop; // positions of coding sequence- start < stop always; strand will tell which is CDS start
  intSiteInfoMap siteData;
  intStrIntMap siteSampleCounts;
  geneData():chromo(),strand(),start(),stop(),siteData(),siteSampleCounts() {}
  geneData(const geneData& o): chromo(o.chromo),strand(o.strand),start(o.start),stop(o.stop),siteData(o.siteData),siteSampleCounts(o.siteSampleCounts) {}
  geneData(std::string c,std::string s,int s1,int s2,const intSiteInfoMap& sd,const intStrIntMap& ssc):chromo(c),strand(s),start(s1),stop(s2),siteData(sd),siteSampleCounts(ssc) {}
};

typedef std::map<std::string, geneData> GeneDataMap;
// we will presume a file with tab-delimited lines that follow this format:
//  0: sampleID
//  1: chromosome
//  2: position
//  3: strand
//  4: tag count
//  5: count of uniq sequence tags that aligned here
//  6: average multiplicity of alignment of tags at this site
//  7: upstream sequence
//  8: downstream sequence
//  9: agCount in the next 10 nt
// 10: aCount in next 10
// 11: gCount in next 10
// 12: site label string
// 13: gene region descriptor
// 14: distance to closest 3'-end
// 15: geneID
// 16: gene start
// 17: gene stop
// 18: gene type
// 19: feature region descriptor
// 20: distance to closest feature 3'-end
// 21: feature ID
// 22: feature start
// 23: feature stop
// 24: feature type

int main(int argc,char *argv[]) {
  if(argc < 2) {
    std::cout << "usage: "<< argv[0] << ": filename\n";
    exit(1);
  }
  GeneDataMap AllGeneData;
  strIntMap AllSamples;
  std::ifstream infile(argv[1],std::ios::in);
  int lineCtr=1,pLineCtr=0;
  std::cerr << "processing input file: " << argv[1] << "..." << std::endl;
  std::ofstream bfile("badAssignedLines.txt",std::ios::out);
  while(!infile.eof()) {
    if(!(lineCtr % 10000)) {
      std::cerr << ".";
      if(!(lineCtr % 2000000)) std::cerr << std::endl;
      else std::cerr.flush();
    }
    char lineIn[1001];
    infile.getline(lineIn,1000);
    strVec f=split()(std::string(lineIn));
    
    if(f.size()==25) {
      ++pLineCtr;
      // get data into readable names and proper type
      std::string geneID(f[15]),chrom(f[1]),str(f[3]),sampID(f[0]);
      int pos(atoi(f[2].c_str())),gStart(atoi(f[16].c_str())),gStop(atoi(f[17].c_str())),tagCount(atoi(f[4].c_str())),agCount[3];
      agCount[0]= atoi(f[9].c_str());
      agCount[1]= atoi(f[10].c_str());
      agCount[2]= atoi(f[11].c_str());
      double avgGcount(atof(f[6].c_str()));
      // first check if the gene has been seen-- if not, then create the record
      if(!AllGeneData.count(geneID)) 
	AllGeneData[geneID]= geneData(chrom,str,gStart,gStop,intSiteInfoMap(),intStrIntMap());
      // next check if this site has been seen; it will be captured simultaneously into both site data sets, so only check one
      if(!AllGeneData[geneID].siteData.count(pos)) {
	AllGeneData[geneID].siteData[pos]= siteInfo(agCount,avgGcount);
	AllGeneData[geneID].siteSampleCounts[pos]= strIntMap();
      }
      // now that we know the record exists, put the count data into the string (sample) to int (count map)
      AllGeneData[geneID].siteSampleCounts[pos][sampID]= tagCount;
      // save the sample names for later
      if(!AllSamples.count(sampID)) AllSamples[sampID]=1;
    } else
      bfile << f.size() << ": " << lineIn << "\n";
    ++lineCtr;
  }
  infile.close();
  bfile.close();
  std::cerr << "done.\n" << pLineCtr << " lines processed out of " << lineCtr << "\n";
  unsigned totalSites=0;
  for(GeneDataMap::iterator it=AllGeneData.begin();it!=AllGeneData.end();++it) totalSites += (*it).second.siteSampleCounts.size();
  std::cerr << AllSamples.size() << " distinct samples\n" << AllGeneData.size() << " distinct genes\n" << totalSites
	    << " distinct sites\n";
  // change of plans-  adding the zeroes made the data structure much larger than was really needed, so instead, just write it out as zeros in the output file
  int geneCtr=1,addCtr=0;
  std::string ofPrefix("joined_site_table");
  if(argc>2) ofPrefix= std::string(argv[2]);
  std::ofstream gfile(std::string(ofPrefix+std::string("_gene.txt")).c_str(),std::ios::out);
  std::ofstream sfile(std::string(ofPrefix+std::string("_site.txt")).c_str(),std::ios::out);
  std::cerr << "Now writing output to files: " << std::string(ofPrefix+std::string("_gene.txt")) << " and "
	    << std::string(ofPrefix+std::string("_site.txt")) << std::endl;
  gfile << "#gene\tchromo\tstart\tstop\tstrand\tsites\n";
  sfile << "#gene\tchromo\tposition\tstrand\tdistToCDSstop\tdistToCDSstart\tavgGcount\tAGcount\tAcount\tGcount";
  for (strIntMap::iterator sIt=AllSamples.begin();sIt!=AllSamples.end();++sIt) sfile << '\t' << (*sIt).first;
  sfile << "\n";
  // put headers into the files  
  for(GeneDataMap::iterator gIt=AllGeneData.begin();gIt!=AllGeneData.end();++gIt) {
    // first write the gene data to the gene file
    gfile << (*gIt).first << '\t' << (*gIt).second.chromo << '\t' << (*gIt).second.start << '\t' << (*gIt).second.stop
	  << '\t' << (*gIt).second.strand << '\t' << (*gIt).second.siteSampleCounts.size() << '\n';
    // now loop through all sites
    intSiteInfoMap::iterator psIt=(*gIt).second.siteData.begin();
    for(intStrIntMap::iterator pIt=(*gIt).second.siteSampleCounts.begin();pIt!=(*gIt).second.siteSampleCounts.end();++pIt) {
      // first write out the stuff about the site itself (position (absolute and 5'- and 3'-relative), multiplicity in genome, agCounts)
      sfile << (*gIt).first << '\t' << (*gIt).second.chromo << '\t' << (*pIt).first << '\t' << (*gIt).second.strand << '\t';
      if((*gIt).second.strand == std::string("+")) sfile << ((*pIt).first - (*gIt).second.stop) << '\t' << ((*pIt).first - (*gIt).second.start);
      else sfile << ((*gIt).second.start-(*pIt).first) << '\t' << ((*gIt).second.stop-(*pIt).first);
      sfile << '\t' << (*psIt).second.gCount << '\t' << (*psIt).second.agCounts[0] << '\t' << (*psIt).second.agCounts[1] << '\t' << (*psIt).second.agCounts[2];
      // now loop through all known samples, writing out either the count, or 0 if the sample isn't present
      for(strIntMap::iterator sIt=AllSamples.begin();sIt!=AllSamples.end();++sIt)
	if(!(*pIt).second.count((*sIt).first)) {
	  //	  (*pIt).second[(*sIt).first]=0;
	  sfile << "\t0";
	  ++addCtr;
	} else
	  sfile << '\t' << (*pIt).second[(*sIt).first];
      sfile << '\n';
      ++psIt;
    }
    if(!(geneCtr % 10)) {
      std::cerr << ".";
      if(!(geneCtr % 2000)) std::cerr << std::endl;
      else std::cerr.flush();
    }
    ++geneCtr;
  }
  sfile.close();
  gfile.close();
  std::cerr << "done.\n" << addCtr << " zero entries added\n";
  return 0;
}
