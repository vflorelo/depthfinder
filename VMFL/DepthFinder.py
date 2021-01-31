#!/opt/intel/intelpython3/bin/python
##!/usr/bin/python3
# Naming conventions:
# Camel case; cum = cumulative ; req = required ; cur = current ; abs = absolute ;
# First global revision:
# Most variables were renamed for consistency
# Summary was re-structured for readability
# Issues persist with methylation percent
# No commits will be performed since the repository is closed
# Victor Flores -> vflorelo@gmail.com 19-11-2019
"""
This script allows an one to estimatate how many fragments can be obtained with a particular pair of enzymes.
Usage:
	./DepthFinder.py PstI MspI 15 50 1000 5 10 Gmax_275_v2_0_no-scaffolds.fasta Soybean 1000 50
"""
import argparse
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *
from Bio.SeqRecord import SeqRecord
import math,datetime,time,copy,sys
import natsort
def analyse(enzymeStart,enzymeEnd,methylationPercent,minBpFragments,maxBpFragments,coverage,missingPercent,genomeFile,speciesName,maxBinValue,binStepValue):
	outfile     = genomeFile[:genomeFile.index(".")]
	dfSummary   = open("Synthesis_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minBpFragments)+"_"+str(maxBpFragments)+".txt",'w')
	dfTable     = open("Table_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minBpFragments)+"_"+str(maxBpFragments)+".txt",'w')
	dfDist      = open("Distrib_Fragments_Length_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(binStepValue)+"_"+str(maxBinValue)+".txt",'w')
	dfTable.write("Sequence name\tLength\tEnzymes\tFragments\tExcluded fragments (<"+str(minBpFragments)+"bp or >"+str(maxBpFragments)+"bp)\tKept fragments\tFragments affected by methylation\tRemaining fragments\tRequired reads (for "+str(coverage)+"X)\n")
	dfDist.write("Sequence\tBin\tNumber of fragments\n")
	restBatch   = RestrictionBatch([enzymeStart])
	if(enzymeEnd != enzymeStart):
		restBatch.add(enzymeEnd)
	numFrags        = 0
	genomeLength    = 0
	cumSizedFrags   = 0
	cumMethFrags    = 0
	cumAbsMethFrags = 0
	cumRealFrags    = 0
	cumReqReads     = 0
	expMutations    = 0
	minCumRealFrags = 0
	minCumReqReads  = 0
	sizeKeptFrags   = 0
	entries         = []
	for record in SeqIO.parse(genomeFile,"fasta"):
		recordLength  = len(record.seq)
		genomeLength += recordLength
		sites    = restBatch.search(record.seq)
		if(enzymeEnd != enzymeStart):
			allCuts   = sites[restBatch.get(enzymeStart)] + sites[restBatch.get(enzymeEnd)]
			cutsStart = set(sites[restBatch.get(enzymeStart)])
			cutsEnd   = set(sites[restBatch.get(enzymeEnd)])
		else:
			allCuts   = sites[restBatch.get(enzymeStart)]
			cutsStart = set(sites[restBatch.get(enzymeStart)])
			cutsEnd   = {}
		allCuts.sort()
		i = 0
		fragList = []
		for values in allCuts[:-1]:
			if len(cutsEnd) > 0:
				if(values in cutsStart and allCuts[i+1] in cutsEnd):
					startPosition = values - 1
					endPosition   = allCuts[i+1] - 1
					fragment      = len(record.seq[startPosition:endPosition])
					if(fragment >= minBpFragments and fragment <= maxBpFragments):
						entry         = SeqRecord(record.seq[startPosition:endPosition],id="Fragment_"+enzymeStart+"-"+enzymeEnd+"_"+str(startPosition)+"-"+str(endPosition)+"_Length="+str(fragment),description=record.description)
						entries.append(entry)
					fragList.append(fragment)
				elif(values in cutsEnd and allCuts[i+1] in cutsStart):
					startPosition = values - 1
					endPosition   = allCuts[i+1] - 1
					fragment      = len(record.seq[startPosition:endPosition])
					if(fragment >= minBpFragments and fragment <= maxBpFragments):
						entry         = SeqRecord(record.seq[startPosition:endPosition],id="Fragment_"+enzymeEnd+"-"+enzymeStart+"_"+str(startPosition)+"-"+str(endPosition)+"_Length="+str(fragment),description=record.description)
						entries.append(entry)
					fragList.append(fragment)
				else:
					pass
			else:
				startPosition = values - 1
				endPosition   = allCuts[i+1] - 1
				fragment      = len(record.seq[startPosition:endPosition])
				if(fragment >= minBpFragments and fragment <= maxBpFragments):
					entry         = SeqRecord(record.seq[startPosition:endPosition],id="Fragment_"+enzymeStart+"-"+enzymeStart+"_"+str(startPosition)+"-"+str(endPosition)+"_Length="+str(fragment),description=record.description)
					entries.append(entry)
				fragList.append(fragment)
			i += 1
		curFrags        = copy.deepcopy(len(fragList))
		curMethFrags    = 0
		curAbsMethFrags = 0
		curRealFrags    = 0
		curReqReads     = 0
		numFrags  += curFrags
		if len(fragList)>0:
			a=1
			b=int(maxBinValue)
			c=int(binStepValue)
			d=int(binStepValue)
			DistFreq={(b+1,max(fragList)):0}
			while c<=b:
				DistFreq[(a,c)]=0
				c+=d
				a+=d
			curSizedFrags = 0
			for lengthValue in fragList:
				if(lengthValue <= minBpFragments or lengthValue >= maxBpFragments):
					curSizedFrags += 1
				for k in DistFreq:
					if (lengthValue >= k[0] and lengthValue <= k[1]):
						DistFreq[k]+=1
					else:
						pass
			for x in natsort.natsorted(DistFreq.items()):
				dfDist.write(record.name+"\t"+str(x[0][0])+"-"+str(x[0][1])+"\t"+str(x[1])+"\n")
			cumSizedFrags   += curSizedFrags
			curMethFrags     = round(((curFrags-curSizedFrags)*methylationPercent)/100)
			cumMethFrags    += curMethFrags
			curAbsMethFrags  = round((curFrags*methylationPercent)/100)
			cumAbsMethFrags += curAbsMethFrags
			curRealFrags     = (curFrags - curSizedFrags) - curMethFrags
			cumRealFrags    += curRealFrags
			curReqReads      = coverage*curRealFrags
			cumReqReads     += curReqReads
			dfTable.write(record.name+"\t"+str(recordLength)+"\t"+enzymeStart+"-"+enzymeEnd+"\t"+str(curFrags)+"\t"+str(curSizedFrags)+"\t"+str(curFrags-curSizedFrags)+"\t"+str(curMethFrags)+"\t"+str(curRealFrags)+"\t"+str(curReqReads)+"\n")
		else:
			dfTable.write(record.name+"\t"+str(recordLength)+"\t"+enzymeStart+"-"+enzymeEnd+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0"+"\t"+"0\n")
	sizeKeptFrags   = numFrags - cumSizedFrags
	expMutations    = round(-0.81+(0.68*(math.log(genomeLength,10))))
	minCumRealFrags = round(cumRealFrags*((100-missingPercent)/100))
	minCumReqReads  = minCumRealFrags * coverage
	SeqIO.write(entries,"Selected_Sequences_%s_%s_%s_%s_%s.fna"%(record.name,enzymeStart,enzymeEnd,str(minBpFragments),str(maxBpFragments)),"fasta")
	dfSummary.write("Species name:\t\t\t\t\t\t\t"+speciesName+"\n")
	dfSummary.write("Sequence name:\t\t\t\t\t\t\t"+genomeFile+"\n")
	dfSummary.write("Sequence length:\t\t\t\t\t\t"+str(genomeLength)+" bp\n")
	dfSummary.write("Enzyme combination:\t\t\t\t\t\t"+enzymeStart+" & "+enzymeEnd+"\n")
	dfSummary.write("Methylation percent:\t\t\t\t\t\t"+str(methylationPercent)+"\n")
	dfSummary.write("Total number of fragments:\t\t\t\t\t"+str(numFrags)+"\n")
	dfSummary.write("Size limits:\t\t\t\t\t\t\t"+str(minBpFragments)+" to "+str(maxBpFragments)+" bp\n")
	dfSummary.write("Fragments outside "+str(minBpFragments)+" and "+str(maxBpFragments)+" bp:\t\t\t\t"+str(cumSizedFrags)+"\n")
	dfSummary.write("Fragments between "+str(minBpFragments)+" and "+str(maxBpFragments)+" bp:\t\t\t\t"+str(sizeKeptFrags)+"\n")
	dfSummary.write("Total fragments affected by methylation:\t\t\t"+str(cumAbsMethFrags)+"\n")
	dfSummary.write("Fragments between "+str(minBpFragments)+" and "+str(maxBpFragments)+" bp affected by methylation:\t"+str(cumMethFrags)+"\n")
	dfSummary.write("Remaining fragments (up to "+str(missingPercent)+"% fragment loss):\t\t\t"+str(minCumRealFrags)+" to "+str(cumRealFrags)+"\n")
	dfSummary.write("Required NGS reads for "+str(coverage)+"X depth:\t\t\t\t"+str(minCumReqReads)+" to "+str(cumReqReads)+" reads\n")
	dfDist.close()
	dfSummary.close()
	dfTable.close()
try:
	enzymeStart        =     sys.argv[1]
	enzymeEnd          =     sys.argv[2]
	methylationPercent = int(sys.argv[3])
	minBpFragments     = int(sys.argv[4])
	maxBpFragments     = int(sys.argv[5])
	coverage           = int(sys.argv[6])
	missingPercent     = int(sys.argv[7])
	genomeFile         =     sys.argv[8]
	speciesName        =     sys.argv[9]
	maxBinValue        = int(sys.argv[10])
	binStepValue       = int(sys.argv[11])
except:
	print("One or more option(s) are missing.")
	print(__doc__)
	sys.exit(1)
analyse(enzymeStart,enzymeEnd,methylationPercent,minBpFragments,maxBpFragments,coverage,missingPercent,genomeFile,speciesName,maxBinValue,binStepValue)
