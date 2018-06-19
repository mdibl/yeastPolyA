#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "/opt/software/internal/yeastPolyA/progs23Nov/util.h"

// presuming that we will get a "*_tgt_aln.txt" file, a tab-delimited text file, in which each line is of the form
//  site<tab>sequence<tab>countOfSeqsAlignedToThisSite<tab>countOfSitesThatThisSeqAlignsTo
// sites have the form chr_position_strand
// sequences have the form
