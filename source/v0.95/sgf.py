# Simple Gappiness Filter to trim alignments (SGF)
#
# This removes sites with too little data to be helpful for phylogenetic analyses.
# 
# It does not distinguish terminal versus internal gaps, and therefore
# it is not intended to assess alignment quality or phylogenetic noise.
# However, used with permissive settings it can remove sites with too little information
# to be worth the computational cost or possibility of artifacts of including low occupancy regions.

# Wayne Maddison

version = "0.95, July 2024"

# example command:
# python sgf.py -input alignment.fas -output trimmedAlignment.fas -siteGT 0.5 -gS -gB -blockSize 5 -blockGT 0.5 -boundary 4
# or
# python3 sgf.py -input alignment.fas -output trimmedAlignment.fas -siteGT 0.5 -gS -gB -blockSize 5 -blockGT 0.5 -boundary 4
# Requires python3

#============================= setting up parameters ===================
inFileName = "infile.fas"
outFileName = "outfile.fas"

#Gappiness filter
filterSiteGappiness = False 
filterBlockGappiness = False

siteGappinessThreshold = 0.5 # A site is considered good (for gappiness) if it is less gappy than this (term or non-term).
minGappyBlockSize = 5 # If in a block of at least this many sites, the first and last site is bad,
blockGappinessThreshold = 0.5 # and the proportion of bad sites is this or above,
minGappyBoundary = 4 # and there are no stretches of this many good sites in a row,

forgiveTaxaWithoutData = True # Forgive taxa without any data (i.e., don't count their gaps).

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
parser.add_argument("-f", help="Forgive taxa without any data (i.e., don't count their gaps).", action='store_true')
parser.add_argument("-not.f", help="Count gaps even in taxa with no data.", dest='f', action='store_false')
parser.set_defaults(f=True)

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
forgiveTaxaWithoutData = args.f


verbose = args.v



if (verbose):
	print("\nSimple Gappiness Filter (SGF) (version ", version, ") with parameters: ")
	if (filterBlockGappiness or filterSiteGappiness):
		print("\n Proportion of gaps for site to be considered gappy: ", args.siteGT)
		if (forgiveTaxaWithoutData):
			print("\n Not counting gaps in taxa with no data.")
		else:
			print("\n Counting gaps even in taxa with no data.")
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
	print("Simple Gappiness Filter (SGF) (v.", version, "): gB=", args.gB, "gS=", args.gS, "siteGT=", args.siteGT, "blockSize=", args.blockSize, "blockGT=", args.blockGT, "boundary=", args.boundary, "f=", args.f)

#============================= Reading the data ===================
#=============================
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


# Checking this is an alignment, i.e. all sequences are the same length
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
	



toDelete = [False for k in range(numChars)]

def anyData(i):
	for k in range(numChars):
		if (sequences[i][k] != "-"):
			return True
	return False

if (forgiveTaxaWithoutData):
	taxonHasData = [False for i in range(numTaxa)]
	numTaxaCounted = 0
	for i in range(numTaxa):
		taxonHasData[i] = anyData(i)
		if (taxonHasData[i]):
			numTaxaCounted += 1
else:
	numTaxaCounted = numTaxa
	
def countTaxon(t):
	if (forgiveTaxaWithoutData):
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

def gappySite(k):
	return siteGappiness[k]>=siteGappinessThreshold

siteGappiness = [0 for k in range(numChars)]
for ic in range(numChars):
	gapCount = 0
	for it in range(numTaxa):
		if (countTaxon(it) and sequences[it][ic] == "-"):
			gapCount+=1
	siteGappiness[ic] = 1.0*gapCount/numTaxaCounted #in case we are to ignore dataless taxa
	if (gapCount == numTaxaCounted): #/if all gaps, delete regardless
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


#============================= Writing sequences without high conflict regions (trimmed) ===================
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
