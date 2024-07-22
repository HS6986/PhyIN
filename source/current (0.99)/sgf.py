# Simple Gappiness Filter to trim alignments (SGF)
# 
# This removes sites with too little data to be helpful for phylogenetic analyses.
# 
# It does not distinguish terminal versus internal gaps, and therefore
# it is not intended to assess alignment quality or phylogenetic noise.
# However, used with permissive settings it can remove sites with too little information
# to be worth the computational cost or possibility of artifacts of including low occupancy regions.
#
# This trimmer was built to be used with PhyIN trimming.

# Wayne Maddison

version = "0.991, July 2024"
# Distributed under an MIT License (see end of file)
# Look for new versions here: https://github.com/wmaddisn/PhyIN

# example command:
# python sgf.py -input alignment.fas -output trimmedAlignment.fas -siteGT 0.5 -gS -gB -blockSize 5 -blockGT 0.5 -boundary 4 -t -1
# or
# python3 sgf.py -input alignment.fas -output trimmedAlignment.fas -siteGT 0.5 -gS -gB -blockSize 5 -blockGT 0.5 -boundary 4 -t -1
#
# Requires python3
# Requires fasta formatted alignment files

#============================= setting up parameters ===================
inFileName = "infile.fas"
outFileName = "outfile.fas"
filterSiteGappiness = False 
filterBlockGappiness = False
siteGappinessThreshold = 0.5 # A site is considered good (for gappiness) if it is less gappy than this (term or non-term).
minGappyBlockSize = 5 # If in a block of at least this many sites, the first and last site is bad,
blockGappinessThreshold = 0.5 # and the proportion of bad sites is this or above,
minGappyBoundary = 4 # and there are no stretches of this many good sites in a row,
taxonCounting = -1 # Counting option for taxa. Options: -1, don't count gaps in taxa without data; 0, count gaps with taxa as is; 1 or higher, count gaps as if this were the total number of taxa

import argparse
parser = argparse.ArgumentParser("Gappy")
parser.add_argument("-input", help="File to be read. Must be DNA/RNA data and in FASTA format, aligned.", default = "infile.fas")
parser.add_argument("-output", help="File to be written.", default = "outfile.fas")
parser.add_argument("-gS", help="Filter by site gappiness.", action="store_true", default = filterSiteGappiness)
parser.add_argument("-gB", help="Filter by block gappiness.", action="store_true", default = filterBlockGappiness)
parser.add_argument("-siteGT", help="Site gappy threshold. Proportion of gaps for site to be considered gappy. Proportion (0 to 1).", type=float, default = siteGappinessThreshold)
parser.add_argument("-blockSize", help="Minimum size of gappy block. Integer.", type=int, default = minGappyBlockSize)
parser.add_argument("-blockGT", help="Block gappy threshold. Proportion of gappy sites for block to be considered gappy. Proportion (0 to 1)", type=float, default = blockGappinessThreshold)
parser.add_argument("-boundary", help="Boundary size. Minimum size of non-gappy block to stop gappy block. Integer.", type=int, default = minGappyBoundary)
tcExplanation = "Taxon counting option (-1, don't count gaps in taxa without data; 0, count gaps with taxa as is; 1 or higher, count gaps as if this were the total number of taxa)."
tcExplanation += " The last option is helpful if processing separate files with different taxa represented, and yet you want them all to apply a gappiness stringency that would be used for the whole set of taxa."
parser.add_argument("-t", help=tcExplanation, type=int, default = minGappyBoundary)

parser.add_argument("-v", help="Verbose.", action="store_true", default = False)
args = parser.parse_args()
inFileName = args.input
outFileName = args.output
filterBlockGappiness = args.gB
filterSiteGappiness = args.gS
siteGappinessThreshold = args.siteGT
minGappyBlockSize = args.blockSize
blockGappinessThreshold = args.blockGT
minGappyBoundary = args.boundary
taxonCounting = args.t
verbose = args.v

if (verbose):
	print("\nSimple Gappiness Filter (SGF) (version ", version, ") with parameters: ")
	if (filterBlockGappiness or filterSiteGappiness):
		print("\n Proportion of gaps for site to be considered gappy: ", args.siteGT)
		if (taxonCounting==-1):
			print("\n Not counting gaps in taxa with no data.")
		elif (taxonCounting==0):
			print("\n Counting gaps even in taxa with no data.")
		elif (taxonCounting>0):
			print("\n Counting gaps as if the total number of taxa were", taxonCounting)
		if (filterSiteGappiness):
			print("\n Filtering individual sites that are too gappy.")
		if (filterBlockGappiness):
			print("\n Filtering blocks of sites that are too gappy:")
			print("   Minimum size of gappy block: ", args.blockSize)
			print("   Proportion of gappy sites for block to be considered gappy: ", args.blockGT)
			print("   Minimum size of non-gappy block to stop gappy block: ",  args.boundary, "\n")
		print("Note: This simple gappiness filter does not distinguish between terminal gaps and internal gaps, nor does it ignore taxa that have no sequence at all.\n")
	else:
		print("\n   NOTHING DONE BY GAPPINESS FILTER: You must indicate whether you want sites (-gS) and/or blocks (-gB) to be filtered")
else:
	if (filterBlockGappiness or filterSiteGappiness):
		print("Simple Gappiness Filter (SGF) (v.", version, "): gB=", args.gB, "gS=", args.gS, "siteGT=", args.siteGT, "blockSize=", args.blockSize, "blockGT=", args.blockGT, "boundary=", args.boundary, "t=", args.t)
	else:
		print("NOTHING DONE BY GAPPINESS FILTER: You must indicate whether you want sites (-gS) and/or blocks (-gB) to be filtered")

#============================= Reading the data ===================
# Read FASTA file
fastaFile = open(inFileName,'r')

# these are the lists that will hold the names and sequences.
names = [] # taxon names
sequences = [] # sequences

countTaxa = 0
currentTaxon = 0
print("Reading file: ", inFileName)
for line in fastaFile: 
	line = line.strip() # strips end-of-line character
	if len(line)>0: #might not be needed, but in case there are blank lines, ignore them
		if line[0] == '>':  #The next taxon!
			names.append(line[1:]) #add this taxon name to the list of taxon names
			sequences.append("") #append an empty string to the list of sequences
			currentTaxon = len(sequences)-1 #remember what taxon number we're now on
			countTaxa +=1
			##print("Reading " + line[1:])
		else: 
			sequences[currentTaxon] += line #not a taxon, therefore sequence, therefore concatenate to the current taxon's sequence
fastaFile.close()
if (verbose):
	print(countTaxa, "sequences in alignment")


# Minimal error checking
if len(sequences) == 0:
	print("ERROR: No sequences read.")
	exit(42)

numChars = len(sequences[0])
for x in sequences:
	if (len(x) != numChars):
		print("ERROR: This appears not to be an alignment; not all taxa have the same sequence length.")
		exit(43)
# perhaps put in here a check that it is DNA data? This is important because of assumption of 4 states + gap		
		
numTaxa = len(names)
if (verbose):
	print("\nNumber of sequences (taxa): ", numTaxa)
	print("Number of sites: ", numChars)
else:
	print("Number of sequences (taxa): ", numTaxa, "; Number of sites: ", numChars)

#..........................
# Setting up information about whether taxa have any data, in case we are forgiving taxa without data (taxonCounting = -1)
def anyData(i):  ## does a taxon have any data at all?
	for k in range(numChars):
		if (sequences[i][k] != "-"):
			return True
	return False

# Gappiness of a column is #gaps/#taxa, but how taxa are counted depends on the taxonCounting option
if (taxonCounting==-1): #not counting gaps in taxa without data; don't even count those taxa
	taxonHasData = [False for i in range(numTaxa)]
	numTaxaCounted = 0
	for i in range(numTaxa):
		taxonHasData[i] = anyData(i)
		if (taxonHasData[i]):
			numTaxaCounted += 1
elif (taxonCounting>0): # counting gaps as if the number of taxa is as specified
	if (taxonCounting<numTaxa):
		print("ERROR: Specified number of taxa in SGF (", taxonCounting, ") is smaller than the number of taxa in the file",numTaxa,")")
		numTaxaCounted = numTaxa
	else:
		numTaxaCounted = taxonCounting
else:
	numTaxaCounted = numTaxa  # counting gaps with the taxa as is
	
def countTaxon(t):
	if (taxonCounting==-1):
		return taxonHasData[t]
	return True


#============================= Gappiness filter ===================
#	CCCCAA--A--A------TTTT--AACCCC
#	CCCCAAG-A--A------TTTT--AACCCC
#	CCCCAAG-A--A------TTTT--AACCCC
#	CCCCAAG----A------------AACCCC
#	       ***********
#
#	In the above, the column after the Gs is the first bad site, and the two following As columns are too narrow to stop the block. However, the Ts
#	columns are 4 in a row, so they stop the block

# This array will record whether the site is marked for deletion
toDelete = [False for k in range(numChars)]

def gappySite(k):
	return siteGappiness[k]>=siteGappinessThreshold

siteGappiness = [0 for k in range(numChars)]
for ic in range(numChars):
	gapCount = 0
	if (taxonCounting>0): #using aspecified number of taxa
		gapCount = taxonCounting - numTaxa  #start with a gap count that is number of taxa not included in file!
		if (gapCount <0):
			gapCount = 0
	for it in range(numTaxa):
		if (countTaxon(it) and sequences[it][ic] == "-"):
			gapCount+=1
	siteGappiness[ic] = 1.0*gapCount/numTaxaCounted #in case we are to ignore dataless taxa
	if (gapCount == numTaxaCounted): #/if all gaps, delete regardless of settings
		toDelete[ic] = True
	if (filterSiteGappiness and gappySite(ic)):
		toDelete[ic] = True		

if (filterBlockGappiness):
	def isGapBlockBoundary(k):
		i = 0
		while (k+i<numChars and i<minGappyBoundary):
			if (gappySite(k+i)):
				return False
			i += 1
		return True

	for blockStart in range(numChars): #let's look for gappy blocks
		if (gappySite(blockStart)): # possible start of gappy block
			boundaryFound = False
			candidateNextBoundary = blockStart+1
			while (candidateNextBoundary<numChars+1 and not boundaryFound): # go ahead until next boundary reached
				if (isGapBlockBoundary(candidateNextBoundary) or candidateNextBoundary == numChars-1):
					boundaryFound = True
					blockEnd = candidateNextBoundary-1

					#blockStart is the potential start of a block; blockEnd is a possible end. If the block is long enough, ask if its blockGappiness is bad
					if (blockEnd-blockStart+1 >= minGappyBlockSize):
						#block is big enough, but is it bad enough?
						badSiteCount = 0
						for k in range(blockStart, blockEnd+1):
							if (gappySite(k)):  # stored as double[] in case criterion shifts, e.g., to average
								badSiteCount +=1
						blockGappiness = 1.0*badSiteCount/(blockEnd-blockStart+1)
						if (blockGappiness >=blockGappinessThreshold):
							for k in range(blockStart, blockEnd+1):
								toDelete[k] = True
				candidateNextBoundary += 1

numDeletedForGappiness = 0
for ic in range(numChars):
	if (toDelete[ic]):
		numDeletedForGappiness += 1
if (verbose):
	print("Number of sites deleted:", numDeletedForGappiness, "; retained:", (numChars-numDeletedForGappiness), "\n")
else:
	print("  Sites deleted:", numDeletedForGappiness, "; retained:", (numChars-numDeletedForGappiness))


#============================= Writing sequences without high gappiness regions (trimmed) ===================
if (verbose):
	print("Writing trimmed sequences to file " + outFileName)
#Sites to delete have been found.  Now to write the output file without them
outputFile = open(outFileName,'w')

for i in range(numTaxa):
	outputFile.write(">")
	outputFile.write(names[i])
	outputFile.write("\n")
	lineLength = 0
	numWritten = 0
	for k in range(numChars):
		if (not toDelete[k]):
			outputFile.write(sequences[i][k])
			numWritten += 1
			lineLength += 1
			if (lineLength % 50 == 0):
				outputFile.write("\n")
		
			
	outputFile.write("\n") 
	
if (verbose):
	print("numWritten ", numWritten)
outputFile.close() 


# MIT License
# 
# Copyright (c) 2024 Wayne Maddison
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.