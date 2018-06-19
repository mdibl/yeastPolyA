#include <iomanip>
#include <iostream>
#include <fstream>

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>

#include "/opt/software/internal/yeastPolyA/progs23Nov/util.h"

class chrFeature {
public:
  std::string group,label;
  int start,stop;
  chrFeature():group(),label(),start(),stop() {}
  chrFeature(int s1,int s2,std::string l,std::string g):group(g),label(l),start(s1),stop(s2) {}
  chrFeature(const chrFeature& o):start(o.start),stop(o.stop),group(o.group),label(o.label) {}
  chrFeature& operator=(const chrFeature& o) { group=o.group; start=o.start; stop=o.stop; label=o.label; return *this; }
  bool operator== (const chrFeature& o) const { return start==o.start && stop==o.stop && label==o.label && group==o.group; }
  bool contains(int pos) const { return pos>=start && pos<=stop; }
  bool operator< (const chrFeature& o) const { return start < o.start; }
  int minDist(int pos) const { if(pos>=start && pos<=stop) return 0; int s1=abs(pos-start); int s2=abs(pos-stop); if(s1<s2) return s1; return s2; }
};

typedef std::vector<chrFeature> cfVector;
typedef std::map<std::string,cfVector> strCfVecMap;
typedef std::map<std::string,strCfVecMap> chrStrVecMap; 

int seqLength=25;
int max5prDist=75,max3prDist=2000,upstreamGeneOverride=300;

int main(int argc,char *argv[]) {
  std::ifstream featFile("/data/internal/Biocore/pcf11_summer/auxiliary/sorted_saccer3_simple_features_with_SUT_CUT_2micron_07Dec2016.txt",std::ios::in);
  //std::ifstream featFile("/Users/jhgraber/data/yeast/saccer3_simple_features.txt",std::ios::in);
  std::cerr << "loading feature file..."; 

  chrStrVecMap FeatureMap;
  while(!featFile.eof()) {
    char line[5001];
    featFile.getline(line,5000);
    strVec f=split()(std::string(line));
    if(f.size()>=6) {
      std::string chr= f[0];
      std::string str= f[4];
      std::string grp= f[1];
      std::string lab= f[5];
      int start= atoi(f[2].c_str());
      int stop= atoi(f[3].c_str());

      if(!FeatureMap.count(chr)) {
	FeatureMap[chr]= strCfVecMap();
	//std::cerr << "Adding chr: " << chr << "\n";
      }
      if(!FeatureMap[chr].count(str)) {
	FeatureMap[chr][str]= cfVector();
	//std::cerr << "Adding chr_str: " << chr << '_' << str << "\n";
      }
      FeatureMap[chr][str].push_back(chrFeature(start,stop,lab,grp));
    }
  }
  featFile.close();
  std::cerr << "Done loading.\n";
  //for(chrStrVecMap::iterator it=FeatureMap.begin();it!=FeatureMap.end();++it)
  //  for(strCfVecMap::iterator it2=(*it).second.begin();it2!=(*it).second.end();++it2)
  //    std::cerr << (*it).first  << ":" << (*it2).first << "\t" << (*it2).second.size() << "\n";

  // now process a "withSeqs siteCount" file
  std::string ifname(argv[1]);
  std::ifstream infile(ifname.c_str(),std::ios::in);
  while(!infile.eof()) {
    char lineIn[5001];
    infile.getline(lineIn,5000);
    strVec f= split()(std::string(lineIn));
    //int ct=0; std::cout << lineIn << "\n"<< ++ct << std::endl;
    if(f.size()>=11)  {
      //UA00e1 chrIII 101129 101179 - 1 CAGTTTCTCTGGGCCACCAACAGCG CAGTTCTGTTTCCCATCTTCATCTT chrIII_101129_101179_-_1 1 1
      //0:sample,1:chr,2:25nt up,3:25nt down,4:strand,5:count,6:upstream25nt,7:downstream25nt,8:label,9:sCount,10:aCount
      std::string chr(f[1]),str(f[4]);
      int site= seqLength+atoi(f[2].c_str());
      int count= atoi(f[5].c_str());
      int scount= atoi(f[9].c_str());
      int acount= atoi(f[10].c_str());

      int dsA=0,dsG=0;
      ///       std::cout << ++ct << "\t" << f[7] << std::endl;
      for(int i=0;i<10;++i) 
	if(f[7][i] == 'A') ++dsA;
 	else if (f[7][i] == 'G') ++dsG;
      std::string dsSeq(f[7]);
      if(FeatureMap.count(chr)) {
 	if(FeatureMap[chr].count(str)) {
 	  int closestFeatureIDX=0;
 	  int closestFeatureDist= FeatureMap[chr][str][0].minDist(site);
 	  int closestGeneIDX= -1;
 	  int closestGeneDist= 100000000,closestGeneMinDist=100000000;
 	  std::string nameString=FeatureMap[chr][str][0].group;
	  if(nameString.length()>4) nameString= nameString.substr(nameString.length()-4,4);
	  if(nameString == std::string("gene")) {
 	    closestGeneIDX= 0;
 	    closestGeneDist= (str==std::string("+"))?abs(site-FeatureMap[chr][str][0].stop):abs(site-FeatureMap[chr][str][0].start);;
 	  }
 	  for(int i=1;i<FeatureMap[chr][str].size();++i) {
	    // currently always taking the gene with the closest 3'-end (dec2016.  this needs to be more subtle
	    //  sites near the 5' end of large genes with upstream tandem genes are being misassigned
	    // gdist is the distance to the closest feature 3'-end (with proper orientation)
	    // fdist is the minimum distance to any part of a feature, and returns 0 if the site is internal to the feature
 	    int gdist= (str==std::string("+"))?abs(site-FeatureMap[chr][str][i].stop):abs(site-FeatureMap[chr][str][i].start);
 	    int fdist= FeatureMap[chr][str][i].minDist(site);
 	    nameString=FeatureMap[chr][str][i].group;
	    if(nameString.length()>4) nameString= nameString.substr(nameString.length()-4,4);
	    if(fdist<closestFeatureDist) {
 	      closestFeatureIDX=i; closestFeatureDist=fdist;
 	    }
	    // changing the critieria assigning to a gene, 22 Dec 2016, as tandem genes are causing problems when only using the closest 3'-end
	    // 
	    //if(nameString==std::string("gene") && (closestGeneIDX == -1 || gdist < closestGeneDist)) {
	    if(nameString==std::string("gene")) {
	      bool newBest=false;
	      if(closestGeneIDX == -1) // no genes encountered yet
		newBest=true;
	      else if(gdist < closestGeneDist && fdist < closestGeneMinDist) // current gene is closer in both 3'-distance and overall distance
		newBest=true;
	      else if(fdist < closestGeneMinDist) { // current gene is closer than best-so-far, but the 3'-end is not closer
		if(fdist == 0) { // internal to the current gene- only reassign if upstream is further than the override distance
		  if(closestGeneDist > upstreamGeneOverride) newBest=true;
		} else { // not internal to the current gene- must be 5', so use the distances
		  if(fdist < max5prDist || closestGeneDist > max3prDist) newBest=true;
		}
	      } else if(gdist < closestGeneDist) { // current gene is closer at 3'-end, but not at the overall distance than current-best
		if(closestGeneMinDist == 0) { // current best is internal to a gene
		  if(gdist <= upstreamGeneOverride) newBest=true; // override if the 3'-end of the upstream gene is "close enough"
		} else { // current best is not internal, but is closer at overall distance than current gene- only a current best 5' can do this
		  if(closestGeneMinDist > max5prDist || gdist > max3prDist) newBest=true;
		}
	      }
	      if(newBest) {
		closestGeneIDX=i; closestGeneDist=gdist; closestGeneMinDist=fdist;
	      }
 	    }
 	  }
 	  //std::cout << chr << "\t" << str << "\t" << site << "\t" << count << "\t" << scount << "\t" << acount
 	  // set a label for the relationship between the site and the identified features
 	  std::string FeatRelPos("f_ext_3"),GeneRelPos("ext_3");
 	  if(FeatureMap[chr][str][closestFeatureIDX].contains(site)) 
 	    FeatRelPos= std::string("f_int");
 	  else if(str==std::string("+") && site<FeatureMap[chr][str][closestFeatureIDX].start)
 	    FeatRelPos= std::string("f_ext_5");
 	  else if(str==std::string("-") && site>FeatureMap[chr][str][closestFeatureIDX].stop)
	    FeatRelPos= std::string("f_ext_5");
	  if(closestGeneIDX < 0)
	    GeneRelPos= std::string("n/a");
	  else if(FeatureMap[chr][str][closestGeneIDX].contains(site)) 
	    GeneRelPos= std::string("int");
	  else if(str==std::string("+") && site<FeatureMap[chr][str][closestGeneIDX].start)
	    GeneRelPos= std::string("ext_5");
	  else if(str==std::string("-") && site>FeatureMap[chr][str][closestGeneIDX].stop)
	    GeneRelPos= std::string("ext_5");
	  // print out results
	  std::cout << f[0] << "\t" << f[1] << "\t" << site << "\t" << f[4] << "\t" << f[5] << "\t" << f[9] << "\t" << f[10] 
		    << "\t" << f[6] << "\t" << f[7] << "\t" << (dsA+dsG) << "\t" << dsA << "\t" << dsG << "\t" << f[8] 
		    << "\t" << GeneRelPos << "\t" << closestGeneDist << "\t";
	  if(closestGeneIDX < 0) 
	    std::cout << "n/a\tn/a\tn/a\tn/a";
	  else
	    std::cout << FeatureMap[chr][str][closestGeneIDX].label << "\t" << FeatureMap[chr][str][closestGeneIDX].start 
		      << "\t" << FeatureMap[chr][str][closestGeneIDX].stop << "\t" << FeatureMap[chr][str][closestGeneIDX].group;
	  std::cout << "\t" << FeatRelPos << "\t" << closestFeatureDist << "\t" << FeatureMap[chr][str][closestFeatureIDX].label 
		    << "\t" << FeatureMap[chr][str][closestFeatureIDX].start << "\t" << FeatureMap[chr][str][closestFeatureIDX].stop 
		    << "\t" << FeatureMap[chr][str][closestFeatureIDX].group << "\n";
	}
      }
    }
  }
  infile.close();
	  
  return 0;
}
