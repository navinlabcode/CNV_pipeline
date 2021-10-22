#!/usr/bin/env python

import sys


def main():

	infilename = sys.argv[1]
	outfilename = sys.argv[2]
	statfilename = sys.argv[3]

	chrominfo = fileToDictionary(sys.argv[4], 0)
	bins = fileToArray(sys.argv[5], 0)
	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	counter = 0
	dups = 0
	totalReads = 0
	prevChrompos = ""
	for x in INFILE:
		arow = x.rstrip().split("\t")
		#thisChrom = arow[2][3:]
		thisChrom = arow[2]
		thisChrom = thisChrom.strip("chr")
		thisChrompos = arow[3]


		if thisChrom.find("_") > -1:
			continue
		if thisChrom.find(".") > -1:
			continue
		if thisChrom.find("*") > -1:
			continue
		if thisChrom.find("M")>-1:
			continue
		if thisChrom.find("s")>-1:
			continue
		if thisChrom == "":
			continue
		if thisChrom == "7d5":
			continue
		if thisChrom == "EBV":
			continue
		if thisChrom == "X":
			thisChrom = "23"
		if thisChrom == "Y":
			thisChrom = "24"

		totalReads += 1
		if thisChrompos == prevChrompos:
			dups += 1
			continue
			
		thisChrominfo = chrominfo[thisChrom]
		thisAbspos = int(thisChrompos) + int(thisChrominfo[2])
		
		counter += 1
		#if counter % 100000 == 0:
		#	print counter
		
		indexUp = len(bins) - 1
		indexDown = 0
		indexMid = int((indexUp - indexDown) / 2.0)

		#print thisChrom, thisChrompos, thisAbspos
		while True:
			#print indexDown, indexMid, indexUp
			if thisAbspos >= int(bins[indexMid][2]):
				indexDown = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexMid
			else:
				indexUp = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexDown

			if indexUp - indexDown < 2:
				break

		#print thisChrom, thisChrompos, thisAbspos, bins[indexDown], bins[indexDown+1]
		binCounts[indexDown] += 1
		prevChrompos = thisChrompos
		
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		OUTFILE.write("\t".join(bins[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")

	binCounts.sort()	
	#print len(binCounts)/2
	bcMedianIndex = len(binCounts)//2
	#Median = 0 
	if isinstance(bcMedianIndex, int):
		bcMedian = binCounts[bcMedianIndex]
	if not isinstance(bcMedianIndex, int):
		bcMedian = (binCounts[int(bcMedianIndex)] + binCounts[int(bcMedianIndex) - 1])//2
	STATFILE.write("TotalReads\tDupsRemoved\tReadsKept\tMedianBinCount\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(dups))
	STATFILE.write("\t")
	STATFILE.write(str(counter))
	STATFILE.write("\t")
	STATFILE.write(str(bcMedian))
	STATFILE.write("\n")
	
	INFILE.close()
	OUTFILE.close()
	STATFILE.close()


def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")

	rd = dict()
#	input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
#		if rd.has_key(id):
		if id in rd:
			#rd[id].append(arow)
			print ("duplicate knowngene id = " + id)
			print ("arow =   " + str(arow))
			print ("rd[id] = " + str(rd[id]))
		else:
			rd[id] = arow
		
	input.close()
	return(rd)


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []

	for i in range(skipFirst):
		input.readline()

	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
		
	input.close()
	return(ra)


if __name__ == "__main__":
	main()
