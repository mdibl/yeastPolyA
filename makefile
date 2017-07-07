makefile:
include configFile.txt  #Including the configuration file allows user to manipulate variables used through makefile

%_alnCount.txt: $(DPATH)%.psl #Creates the _alnCount.txt file from the .psl file
	echo $(DPATH) $<
	grep seq $< | perl -ane "if(@F[0]>=(@F[10]*$(FAC)) && @F[17]<2) { print \"@F[9]\n\";} " | uniq -c | perl -ane "print \"@F[1]\t@F[0]\n\"; " | sort > $(outDir)$@

%_tgtCount.txt: $(DPATH)%.psl #Creates the _tgtCount.txt  file from the .psl file
	echo $(DPATH) $<
	grep seq $< | perl -ane "if(@F[0]>=(@F[10]*$(FAC)) && @F[17]<2) { if(@F[8] eq \"-\") { print \"@F[13]_\"; print (@F[16]+@F[11]); print \"_-\n\";} else { print \"@F[13]_\"; print (@F[15]-@F[11]); print \"_+\n\";} } " | sort | uniq -c | perl -ane "print \"@F[1]\t@F[0]\n\"; " > $(outDir)$@ 

%_tgt_aln.txt: $(DPATH)%.psl $(outDir)%_alnCount.txt $(outDir)%_tgtCount.txt #Creates the _tgt_aln.txt file from the _alnCount.txt and _tgtCount.txt files
	echo $(DPATH) $<
	grep seq $< | perl -ane "if(@F[0]>=(@F[10]*$(FAC)) && @F[17]<2) { if(@F[8] eq \"-\") { print \"@F[9]\t@F[13]_\"; print (@F[16]+@F[11]); print \"_-\n\";} else { print \"@F[9]\t@F[13]_\"; print (@F[15]-@F[11]); print \"_+\n\";} } " | sort -k 2 | join -1 2 - $*_tgtCount.txt | perl -ane "print \"@F[1]\t@F[0]\t@F[2]\n\"; " | sort | join - $*_alnCount.txt | perl -ane "print \"@F[1]\t@F[0]\t@F[2]\t@F[3]\n\"" > $(outDir)$@

%_siteCount.txt: $(outDir)%_tgt_aln.txt $(CPATH)CondenseTagsToSiteCount #Creates the _siteCounts.txt file from the CondenseTagsToSiteCount code with _tgt_aln.txt as an arguement
	$(CPATH)CondenseTagsToSiteCount $< > $(outDir)$@

%_siteCount.bed: $(outDir)%_siteCount.txt #Converts the _siteCount.txt file to _siteCount.bed via commandline command
	perl -ne "s/_\+/_#/; s/_\-/_+/; s/_#/_-/; s/_/\t/g; print; " $< | perl -ane "print \"@F[0]\t\"; print (@F[1]-$(PAD)); print \"\t\"; print (@F[1]+$(PAD)); print \"\tsamplename\t@F[3]\t@F[2]\t@F[4]\t@F[5]\n\"; " > $(outDir)$@

%_siteCount_withSeqs.txt: $(outDir)%_siteCount.bed $(auxDir)$(faFile) $(CPATH)faGetBedSequences $(auxDir)$(faiFile) #Creates the _siteCount_withSeqs.txt file by running the faGetBedSequences code with faFile and _siteCount.bed as arguements, that is then piped to a perl command the modifies the datafile and then produces the final output file
	$(CPATH)faGetBedSequences $(auxDir)$(faFile) $< | perl -ne "s/_/\t/g; print" | perl -ane "print \"@F[0]\t@F[1]\t@F[2]\t@F[3]\t@F[4]\t@F[5]\t\"; print substr(@F[6],0,$(PAD)); print \"\t\"; print substr(@F[6],$(PAD),2*$(PAD)); print \"\t@F[7]_\"; print ((@F[2]+@F[3])/2); print \"_@F[4]_@F[5]\t@F[12]\t@F[13]\n\"; " > $(outDir)$@

%_assigned.txt: $(outDir)%_siteCount_withSeqs.txt $(auxDir)$(sortedFeatures) $(CPATH)AssignSitesToGenes #Creates the _assigned.txt file by running the AssignSitesToGenes code with _withSeqs.txt and sortedFeatures as arguements
	$(CPATH)AssignSitesToGenes $< $(auxDir)$(sortedFeatures) > $(outDir)$@  #Need to change the directory that the sortedFeatures file is in, in the c++ code as well

CollectSiteCountsOutputFiles: $(outDir)%_assigned.txt $(CPATH)CollectSiteCounts #Creates the two output files from the CollectSiteCounts code with _assigned.txt as the cole arguement
	$(CPATH)CollectSiteCounts $< > $(outDir)$@



