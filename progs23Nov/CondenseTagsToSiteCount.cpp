//
//  CondenseTagsToSiteCount.cpp
//  
//
//  Created by Joel Graber on 3/24/16.
//
// This program will take the post-processed alignment of polyA site tags and condense all tags
//    that hit the same site to one count, while also recording the, number of tags that hit the site
//    and the maximum ambiguity of the site, as defined by the number of alternative alignments of each
//    tag that aligns.

// The input file is presumed to have lines of the following form:
// chrIII_100708_+	seq_1107896_1	1	1
// chrIII_101014_+	seq_154399_2	1	1
// ...
// Column 1 identifies the site, by chromosome ID (string), position (integer), and strand (+/-).
//    NOTE: tags are presumed to be on the antisense from the originating transcript
// Column 2 has tag information: seq_ + rank (integer)  + _ + count, where the presumption is that
//    prior to alignment, tags were reduced to unique sequences, counted and ranked, with the most
//    abundant tag at rank 1.  The second number (count) is what will be aggregated in the count
// Column 3 is the number of tags that hit this specific location (including orientation)
// Column 4 is tne number of accepted alignments for this sequence tag

// It is presumed that the alignments were prefiltered for low-quality matches before running this program
//    so all tags and counts will be used.
//

// 14 Dec 2016; modifying the program to calculate an average ambiguous count for each site rather than just
//    accepting the maximum as before.  This was found to be too conservative and pessimistic.
//    We will use the counts (third part of the sequence label) as a weighting term on the ambiguous counts for
//    each sequence that aligns to the given site.
//    The principle modification is that the last element of the count vector in the strIntMap will now
//    hold the weighted sum of ambigous counts, so it will need to be divided by the total counts (which is
//    the sum of the weights) on output to the new file.

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include "/Users/knaggert/Desktop/pcf11_summer/progs23Nov/util.h"

// use STL maps to hold the pertinent data
typedef std::vector<int> iVec;
typedef std::map<std::string,iVec,std::less<std::string> > strIntMap;

int main(int argc,char *argv[]) {
    if(argc<2) {
        std::cout << "usage: " << argv[0] << " inputfile\n";
        exit(0);
    }
    strIntMap siteCount;
    std::ifstream infile(argv[1],std::ios::in);
    int nCount=0,aCount=0;

    while (!infile.eof()) {
        char lineIn[1001];
        infile.getline(lineIn,1000);
        strVec f=split()(std::string(lineIn));
        if(f.size()==4) {
            std::string siteID= f[0];
            int tagCount= atoi(f[2].c_str());
            int ambCount= atoi(f[3].c_str());
            int pCount= f[1].find('_',4)+1;
            int count= atoi(f[1].substr(pCount,f[1].length()-pCount).c_str());
            //std::cout << siteID << '\t' << f[1] << '\t' << count << '\t' << tagCount << '\t' << ambCount << '\n';
            if(!siteCount.count(siteID)){
                ++nCount;
                //std::cout << "Adding record for " << siteID << '\t' << nCount << '\t' << aCount << std::endl;
                siteCount[siteID]=iVec(3,0); 
                siteCount[siteID][0]= count;
                siteCount[siteID][1]= tagCount;
                siteCount[siteID][2]= count*ambCount;
            } else {
                ++aCount;
                siteCount[siteID][0]+= count;
		siteCount[siteID][2]+= count*ambCount;
                //if(ambCount>siteCount[siteID][2])
                //    siteCount[siteID][2]= ambCount;
            }
        }
    }
    infile.close();
    
    for(strIntMap::iterator it=siteCount.begin();it!=siteCount.end();++it) {
        std::cout << (*it).first << '\t' << (*it).second[0] << '\t' << (*it).second[1] << '\t'
		  << (static_cast<double>((*it).second[2])/static_cast<double>((*it).second[0])) << '\n';
    }
    
    return 0;
}
