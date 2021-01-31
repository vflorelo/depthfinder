#!/usr/bin/env python3

"""
This script allows an one to estimatate how many fragments can be obtained with a particular pair of enzyme.  

Usage:
	./DepthFinder.py PstI MspI 15 50 1000 5 10 Gmax_275_v2_0_no-scaffolds.fasta Soybean 1000 50

"""

import argparse
#parser = argparse.ArgumentParser()

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *
from Bio.SeqRecord import SeqRecord
import math,datetime,time,copy,sys
import natsort

def analyse(enzymeStart,enzymeEnd,meth,minbp,maxbp,cov,mperc,genome,species,max_b_bin,step_d_bin):
	#print("Current date and time: " + str(datetime.datetime.now()) + "\n")

	methylationPercent=int(meth)
	minBpFragments=int(minbp)
	maxBpFragments=int(maxbp)
	coverage=int(cov)
	missingPercent=int(mperc)

	dot = genome.index(".")
	ext = genome[dot+1:]
	outfile = genome[:dot]
	
	#Creation of output files:
	resume = open("Synthesis_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minbp)+"_"+str(maxbp)+".txt", 'w')
	table = open("Table_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+str(minbp)+"_"+str(maxbp)+".txt", 'w')
	distrib = open("Distrib_Fragments_Length_GBS_"+outfile+"_"+enzymeStart+"_"+enzymeEnd+"_"+step_d_bin+"_"+max_b_bin+".txt", 'w')

	table.write("Species Name\tSequence Name\tSequence Name\tSeq Length\tEnzyme combination\tNb cutting sites\tNb fragments\t%Methylation\tNb fragments affected by methylation\tRemaining fragments\tNb fragments affected by site mutation\tRemaining fragments\tSize limits\tNb fragments discarded\tRemaining fragments\tDepth of coverage (X)\tRequired read counts\tMissing fragments\tRequired read counts\n")
	distrib.write("Sequence\tBin\tNumber of fragments\n")
				
	#Loading restriction enzymes
	rb = RestrictionBatch([enzymeStart])
	if(enzymeEnd != enzymeStart):
		rb.add(enzymeEnd)
	#print("rb: ",rb)
		
	tot_len=0
	nb_cuts=0
	nb_frag=0
	nb_methyl=0
	rem_mfrag_meth=0
	tot_formula=0
	rem_frag_formula=0
	nb_mut=0
	nb_sizelim=0
	toRemAll=0
	smallCutsAll=0
	frag_cov=0
	missFragAll=0
	remFragMiss=0
	
	#Analyze each sequence (chromosomes, scaffolds) in the genome file
	entries=[]
	for record in SeqIO.parse(genome, "fasta"):
		#print(record.name)
		
		tot_len+=len(record.seq)

		#Searching restriction sites
		sites = rb.search(record.seq)
		#print("Dictionnaire 'sites' contenant les positions genomiques du/des site(s) de restriction:")
		#print(sites)
		#for x in sites:
			#print(x,len(sites[x]))

		cuts=0;
		for key, value in sites.items():
			#key: each restriction enzyme
			#value: list of restriction site position in the sequence
			#print(key,value)
			cuts += len(value)
		#print("Number of all cutting sites found : " + str(cuts))
			
		#In the case of a 2 enzymes analysis, we merge all the positions that those 2 enzymes cuts
		#We use Python sets which are more faster for searching a value in
		if(enzymeEnd != enzymeStart):
			allCuts = sites[rb.get(enzymeStart)] + sites[rb.get(enzymeEnd)]
			cutsStart = set(sites[rb.get(enzymeStart)])
			cutsEnd = set(sites[rb.get(enzymeEnd)])
		else:
			allCuts = sites[rb.get(enzymeStart)]
			cutsStart = set(sites[rb.get(enzymeStart)])
			cutsEnd = {}
		
		#Put the positions in order
		allCuts.sort()
		#print("Liste de toutes les positions genomiques de sites de restriction (valeur de allCuts): ",allCuts)
		#print("Nombre de valeurs dans allcuts: ",len(allCuts))
		
		#print("cutsStart: ",cutsStart,len(cutsStart))
		#print("cutsEnd: ",cutsEnd,len(cutsEnd))
			
		#print("\nFrom the list of all cutting sites, scan for site pair where the first site belong to startEnzyme and the second site belong to endEnzyme\n")
		#print("Liste de toutes les positions genomiques de sites de restriction (valeur de allCuts): ",allCuts)

		i = 0
		#entries=[]
		listlength = []
		for values in allCuts[:-1]:# We don't process the last cut
			#print(values)
			if len(cutsEnd) > 0:
				if(values in cutsStart and allCuts[i+1] in cutsEnd):# or (values in cutsEnd and allCuts[i+1] in cutsStart):
					beginingPosition = values - 1
					endPosition = allCuts[i+1] - 1
					fragment=len(record.seq[beginingPosition:endPosition])
					entry=SeqRecord(record.seq[beginingPosition:endPosition],id="Fragment_"+enzymeStart+"-"+enzymeEnd+"_"+str(beginingPosition)+"-"+str(endPosition)+"_Length="+str(fragment),description=record.description)
					entries.append(entry)
					listlength.append(fragment)
					#print(beginingPosition,endPosition,fragment)
				elif(values in cutsEnd and allCuts[i+1] in cutsStart):
					beginingPosition = values - 1
					endPosition = allCuts[i+1] - 1
					fragment=len(record.seq[beginingPosition:endPosition])
					entry=SeqRecord(record.seq[beginingPosition:endPosition],id="Fragment_"+enzymeEnd+"-"+enzymeStart+"_"+str(beginingPosition)+"-"+str(endPosition)+"_Length="+str(fragment),description=record.description)
					entries.append(entry)
					listlength.append(fragment)
					#print(beginingPosition,endPosition,fragment)
				else:
					pass
			else:
				beginingPosition = values - 1
				endPosition = allCuts[i+1] - 1
				fragment=len(record.seq[beginingPosition:endPosition])
				entry=SeqRecord(record.seq[beginingPosition:endPosition],id="Fragment_"+enzymeStart+"-"+enzymeStart+"_"+str(beginingPosition)+"-"+str(endPosition)+"_Length="+str(fragment),description=record.description)
				entries.append(entry)
				listlength.append(fragment)
				#print(beginingPosition,endPosition,fragment)
				
			i += 1
		#print("Nombre de fragments total traites: ",i)


		#The length of listlength give the number of good cutting sites and fragments
		#In the case of one enzyme, this is all cutting sites found
		#In the case of two enzymes, we keep only paired cutting sites and corresponding fragment obtained
		#print("Nombre de fragments avec les bons sites de restriction a chaque extremite: ",len(listlength),"\n")

		#Sum of cuts over all sequences
		gcuts=copy.deepcopy(len(listlength))
		nb_cuts+=gcuts

		#Number of fragments found. Normally, if we cut a rope 3 times, we will have 4 pieces. But here, there are no cutting sites
		#at the beginning nor at the end so we remove the first and the last fragments (they won't be processed).
		#That's because if we have 3 restriction sites, we will process 2 fragments.
		frag=copy.deepcopy(len(listlength))
		nb_frag+=frag

		#File containing fragments
		#SeqIO.write(entries,"Selected_Sequences_%s_%s_%s_%s_%s.fna"%(record.name,enzymeStart,enzymeEnd,str(minbp),str(maxbp)),"fasta")

		#Section for creation of a distribution table for the frequency of fragments based on their length
		if len(listlength)>0:
			a,b,c,d=1,int(max_b_bin),int(step_d_bin),int(step_d_bin)
			# 'a' is the 1st value of the bin
			# 'b' the max length before last bin,
			# 'c' is the 2nd value of the bin
			# 'd' walking value
			# 'c' should be equal to 'd'.

			# here we entre the last bin based on the lenght above max_b_bin
			DistFreq={(b+1,max(listlength)):0}
			while c<=b:
				DistFreq[(a,c)]=0
				c+=d
				a+=d

			for values in listlength:
				for k in DistFreq:
					if (values >= k[0] and values <= k[1]):
						DistFreq[k]+=1
					else:
						pass

			#print("Voici DistFreq:")
			#print(natsort.natsorted(DistFreq.items()))
			
			#Writing the distribution length file
			for x in natsort.natsorted(DistFreq.items()):
				distrib.write(record.name+"\t"+str(x[0][0])+"-"+str(x[0][1])+"\t"+str(x[1])+"\n")

			#Calculate the number of fragment according to methylation
			cutsMeth = frag - ((frag * methylationPercent) / 100)
			nb_methyl+=round(((frag * methylationPercent) / 100))
			rem_mfrag_meth+=round(cutsMeth)

			#Calculate the number of fragments according to site mutation
			formula = round(-0.81 + 0.68 * (math.log(len(record.seq), 10)))
			cutsFormula = cutsMeth - formula
			tot_formula+=round(formula)
			rem_frag_formula+=round(cutsFormula)

			#Calculate the number of fragments of the right length.
			toRemove = 0
			for values in listlength:
				if(values <= minBpFragments or values >= maxBpFragments):
					toRemove += 1
			smallCuts = cutsFormula - toRemove;

			#Fragments to remove
			toRemAll+=toRemove
			#Remianing fragments
			smallCutsAll+=round(smallCuts)

			#Calculate the number of frgaments according to depth of coverage
			cutsCoverage = round(smallCuts) * coverage
			frag_cov+=cutsCoverage
	
			#Calculate the number of fragments according to the percentage of missing sites (missingPercent)
			missFrag=round(((cutsCoverage * missingPercent) / 100))
			cutsMissing = cutsCoverage - missFrag
			missFragAll+=missFrag
			remFragMiss+=cutsMissing
				
			table.write(species + "\t" + genome + "\t" + record.name + "\t" + str(len(record.seq)) + "\t" + enzymeStart + "-" + enzymeEnd + "\t" + str(gcuts) + "\t" + str(frag) + "\t" + str(methylationPercent) + "\t" + str(round((frag * methylationPercent)/100)) + "\t" + str(int(round(cutsMeth))) + "\t" + str(round(formula)) + "\t" + str(int(round(cutsFormula))) + "\t" + str(minBpFragments)+"_" + str(maxBpFragments) + "\t"+ str(toRemove) + "\t" +  str(int(round(smallCuts))) + "\t" + str(coverage) + "\t" + str(int(round(cutsCoverage))) + "\t" + str(round(((cutsCoverage * missingPercent) / 100))) + "\t" + str(int(round(cutsMissing))) + "\n")

		else:
			table.write(species + "\t" + genome + "\t" + record.name + "\t" + str(len(record.seq)) + "\t" + enzymeStart + "-" + enzymeEnd + "\t" + "0" + "\t" + "0" + "\t" + str(methylationPercent) + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\t" + str(minBpFragments)+"_" + str(maxBpFragments) + "\t"+ "0" + "\t" +  "0" + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\n")

	#File containing fragments
	SeqIO.write(entries,"Selected_Sequences_%s_%s_%s_%s_%s.fna"%(record.name,enzymeStart,enzymeEnd,str(minbp),str(maxbp)),"fasta")

	resume.write("Species Name\t:\t"+species+"\n")
	resume.write("Sequence Name\t:"+genome+"\n")
	resume.write("Seq\t:\tWhole Genome"+"\n")
	resume.write("Length\t:\t"+str(tot_len)+"\n")
	resume.write("Enzyme combination\t:\t"+enzymeStart + "-" + enzymeEnd+"\n")
	resume.write("Nb cutting sites\t:\t"+str(nb_cuts)+"\n")
	resume.write("Nb fragments\t:\t"+str(nb_frag)+"\n")
	resume.write("% Methylation\t:\t"+str(meth)+"\n")
	resume.write("Nb fragments affected by methylation\t:\t"+str(int(nb_methyl))+"\n")
	resume.write("Remaining fragments\t:\t"+str(int(rem_mfrag_meth))+"\n")
	resume.write("Site mutation\t:\t"+str(tot_formula)+"\n")
	resume.write("Remaining fragments\t:\t"+str(rem_frag_formula)+"\n")
	resume.write("Size limits\t:\t"+str(minBpFragments)+"-" + str(maxBpFragments)+"\n")
	resume.write("Nb fragments discarded\t:\t"+str(toRemAll)+"\n")
	resume.write("Remaining fragments\t:\t"+str(smallCutsAll)+"\n")
	resume.write("Depth of coverage (X)\t:\t"+str(cov)+"\n")
	resume.write("Required read counts\t:\t"+str(frag_cov)+"\n")
	resume.write("Missing fragments\t:\t"+str(missFragAll)+"\n")
	resume.write("Required read counts\t:\t"+str(remFragMiss)+"\n")

	distrib.close()
	resume.close()
	table.close()

try:
	enzymeStart = sys.argv[1]
	enzymeEnd = sys.argv[2]
	meth = int(sys.argv[3])
	minbp = int(sys.argv[4])
	maxbp = int(sys.argv[5])
	cov = int(sys.argv[6])
	mperc = int(sys.argv[7])
	genome = sys.argv[8]
	species = sys.argv[9]
	max_b_bin = sys.argv[10]
	step_d_bin = sys.argv[11]
except:
	print("One or more option(s) are missing.")
	print(__doc__)
	sys.exit(1)

analyse(enzymeStart,enzymeEnd,meth,minbp,maxbp,cov,mperc,genome,species,max_b_bin,step_d_bin)
