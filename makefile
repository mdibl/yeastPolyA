# Script Name: makefile
# Author: Kristoph S. N. Naggert, Thursday August 10th 2017
# Type: Makefile

# Arguments: Name of Program You Want to Compile or Name of Object File You Want to Create

# Execution: make (ObjectFile or Program)

# Overview:  The makefile builds all of the object files and executable programs in the pcf11 project workflow.

# Variables: No variables.

include /opt/software/internal/yeastPolyA/configFile.txt

%.o: %.cpp #Creates the .o file for any given .cpp file
	g++ -O3 -c $<

CondenseTagsToSiteCount: CondenseTagsToSiteCount.o util.o logFile.o #Compiles the CondenseTagsToSiteCount code
	g++ -O3 -o $@ $^ -lgsl -lgslcblas

faGetBedSequences: faGetBedSequences.o util.o logFile.o basicSeq.o fastaSeq.o filter.o #Compiles the faGetBedSequences code
	g++ -O3 -o $@ $^ -lgsl -lgslcblas

AssignSitesToGenes: AssignSitesToGenes.o util.o logFile.o #Compiles the AssignSitesToGenes code
	g++ -O3 -o $@ $^ -lgsl -lgslcblas

CollectSiteCounts: collectSiteCounts.o util.o logFile.o #Compiles the CollectSiteCounts code
	g++ -O3 -o $@ $^ -lgsl -lgslcblas

processAllPolyASites: processAllPolyASites.o util.o logFile.o #Compiles the processAllPolyASites code
	g++ -O3 -o $@ $^ -lgsl -lgslcblas
